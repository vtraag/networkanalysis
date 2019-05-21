package cwts.networkanalysis;

import java.util.Random;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;

/**
 * Fast labeled target cover local moving algorithm.
 *
 * <p>
 * The algorithm moves nodes from one cluster to another cluster to improve the
 * quality function as much as possible. The algorithm starts with a queue of
 * all nodes. It removes the node from the front queue, and moves it to the
 * cluster that maximizes the quality function. The node remains in the same
 * cluster if no alternative cluster is better. If the cluster has changed, all
 * its neighbors from the old cluster that are not yet in the queue are
 * appended to the end of the queue. The algorithm continues this process until
 * the queue is empty.
 * </p>
 *
 * <p>
 * The initial order of the nodes in the queue is randomized.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 */
public class FastLocalMovingAlgorithmLabeledTargetCover extends IterativeCPMClusteringAlgorithmLabeledTargetCover
{
    protected Random random;
    
    /**
     * Constructs fast labeled target cover local moving algorithm.
     *
     * <p>
     * Uses default resolution. Uses a random number generator with a random
     * seed.
     * </p>
     */
    public FastLocalMovingAlgorithmLabeledTargetCover()
    {
        this(DEFAULT_RESOLUTION, new Random());
    }

    /**
     * Constructs fast labeled target cover local moving algorithm using specified
     * random number generator.
     *
     * <p>
     * Uses default resolution.
     * </p>
     *
     * @param random Random number generator
     */
    public FastLocalMovingAlgorithmLabeledTargetCover(Random random)
    {
        this(DEFAULT_RESOLUTION, random);
    }

    /**
     * Constructs fast labeled target cover local moving algorithm for specified
     * resolution parameter and specified random number generator.
     *
     * @param resolution Resolution parameter
     * @param random Random number generator
     */
    public FastLocalMovingAlgorithmLabeledTargetCover(double resolution, Random random)
    {
        super(resolution);
        this.random = random;
    }

    /**
     * Constructs fast labeled target cover local moving algorithm for specified
     * resolution parameter and specified random number generator.
     *
     * @param resolution Resolution parameter
     * @param random Random number generator
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public FastLocalMovingAlgorithmLabeledTargetCover(double resolution, Random random, Cover targetCover, double targetWeight)
    {
        super(resolution, targetCover, targetWeight);
        this.random = random;
    }
    
    /**
     * Constructs fast labeled target cover local moving algorithm for specified
     * resolution parameter and specified random number generator.
     *
     * @param resolution Resolution parameter
     * @param nIterations Number of iterations
     * @param random Random number generator
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public FastLocalMovingAlgorithmLabeledTargetCover(double resolution, int nIterations, Random random, Cover targetCover, double targetWeight)
    {
        super(resolution, nIterations, targetCover, targetWeight);
        this.random = random;
    }

    /**
     * Gets random number generator.
     *
     * @return Random number generator
     */
    public Random getRandom()
    {
        return random;
    }

    /**
     * Sets random number generator.
     *
     * @param random Random number generator
     */
    public void setRandom(Random random)
    {
        this.random = random;
    }

    /**
     * Runs the fast labeled target cover local moving algorithm.
     *
     * <p>
     * The algorithm moves nodes from one cluster to another cluster to improve
     * the quality function as much as possible. The algorithm starts with a
     * queue of all nodes. It removes the node from the front queue, and moves
     * it to the cluster that maximizes the quality function. The node remains
     * in the same cluster if no alternative cluster is better. If the cluster
     * has changed, all its neighbors from the old cluster that are not yet in
     * the queue are appended to the end of the queue. The algorithm continues
     * this process until the queue is empty.
     * </p>
     *
     * <p>
     * It may be beneficial to move a node to a cover to which it is assigned in
     * the target cover. Because of this, the number of required clusters may
     * actually be larger than the number of nodes, since the cover can be an
     * arbitrary number.
     * </p>
     *
     * <p>
     * No empty clusters are ever removed. This would result in relabeling of
     * the clusters, which may decrease the agreement with the target cover.
     * </p>
     *
     * @return Boolean indicating whether the clustering has changed
     */
    protected boolean improveClusteringOneIteration(Network network, Clustering clustering)
    {
        boolean update;
        boolean[] stableNode, addedClusterForConsideration;
        double maxQualityValue, qualityValue;
        double[] clusterWeight, edgeWeightPerCluster, qualityValuePerConsideredCluster;
        int bestCluster, currentCluster, i, j, k, l, nClustersForConsideration, nUnstableNodes;
        int[] clusterForConsideration, nNodesPerCluster, nodeOrder;
        int nextUnusedCluster;
        Cover nodesPerCover = this.targetCover.getNodesPerCover();
        
        update = false;

        /*
         * Possibly a node should be assigned to one of its covers. Since the
         * cover is an arbitrary number, we need to ensure to allocate
         * sufficient possible clusters.
         */
        int maxNClusters = Math.max(clustering.nClusters, targetCover.nCovers) + network.nNodes;

        // Calculate total node weight and number of nodes per cluster.
        clusterWeight = new double[maxNClusters];
        nNodesPerCluster = new int[maxNClusters];
        for (i = 0; i < network.nNodes; i++)
        {
            clusterWeight[clustering.clusters[i]] += network.nodeWeights[i];
            nNodesPerCluster[clustering.clusters[i]]++;
        }

        /*
         * Determine the next cluster that is unused. We do this in a slightly
         * different way than ordinary, because of the labeled target cover.
         * Normally, we could take any empty cluster. However, we now want to
         * take an empty cluster that is not in used by any cover, so that other
         * nodes from that cover may be assigned to that cluster.
         */
        nextUnusedCluster = maxNClusters;
        for (i = 0; i < maxNClusters; i++)
            // Pick the first unused cluster that is not covered 
            if (nNodesPerCluster[i] == 0 && ((i < nodesPerCover.nNodes && nodesPerCover.getNCovers(i) == 0) || i >= nodesPerCover.nNodes))
            {
                nextUnusedCluster = i;
                break;
            }

        // Generate random permutation of nodes.
        nodeOrder = cwts.util.Arrays.generateRandomPermutation(network.nNodes, random);

        edgeWeightPerCluster = new double[maxNClusters];
        clusterForConsideration = new int[maxNClusters];
        stableNode = new boolean[network.nNodes];
        nUnstableNodes = network.nNodes;
        addedClusterForConsideration = new boolean[maxNClusters];
        qualityValuePerConsideredCluster = new double[maxNClusters];
        i = 0;
        double improvement = 0.0;
        do
        {
            j = nodeOrder[i];
            currentCluster = clustering.clusters[j];

            // Remove the currently selected node from its current cluster.
            clusterWeight[currentCluster] -= network.nodeWeights[j];

            /*
             * Decrease the number of node of the cluster, and check for the
             * next empty clusters.
             */
            nNodesPerCluster[currentCluster]--;
            if (nNodesPerCluster[currentCluster] == 0)
            {
                // If the currentCluster is smaller and not used, set that as next cluster
                if (currentCluster < nextUnusedCluster && ((currentCluster < nodesPerCover.nNodes && nodesPerCover.getNCovers(currentCluster) == 0) || currentCluster >= nodesPerCover.nNodes))
                    nextUnusedCluster = currentCluster;
            }
            

            // Add the empty cluster for consideration
            clusterForConsideration[0] = nextUnusedCluster;
            addedClusterForConsideration[clusterForConsideration[0]] = true;
            nClustersForConsideration = 1;
            
            /*
             * Identify the neighboring clusters of the currently selected node,
             * that is, the clusters with which the currently selected node is
             * connected.
             */
            for (k = network.firstNeighborIndices[j]; k < network.firstNeighborIndices[j + 1]; k++)
            {
                l = clustering.clusters[network.neighbors[k]];
                if (!addedClusterForConsideration[l])
                {
                    clusterForConsideration[nClustersForConsideration] = l;
                    addedClusterForConsideration[l] = true;
                    nClustersForConsideration++;
                }
                edgeWeightPerCluster[l] += network.edgeWeights[k];
            }

            // Also add covers of node as clusters for consideration
            for(Int2DoubleMap.Entry entry : targetCover.getCoverIterator(j))
            {
                l = entry.getIntKey();
                if (!addedClusterForConsideration[l])
                {
                    clusterForConsideration[nClustersForConsideration] = l;
                    addedClusterForConsideration[l] = true;
                    nClustersForConsideration++;
                }
            }
            
            /*
             * For each cluster for consideration, calculate the increment of
             * the quality function obtained by moving the currently selected
             * node to the neighboring cluster. Determine the cluster for which
             * the increment of the quality function is largest. The currently
             * selected node will be moved to this optimal cluster. In order to
             * guarantee convergence of the algorithm, if the old cluster of the
             * currently selected node is optimal but there are also other
             * optimal clusters, the currently selected node will be moved back
             * to its old cluster.
             */
            bestCluster = currentCluster;
            double currentQualityValue = edgeWeightPerCluster[currentCluster] - network.nodeWeights[j] * clusterWeight[currentCluster] * resolution
                              + targetCover.getCoverWeight(j, currentCluster) * targetWeight;
            maxQualityValue = currentQualityValue;
            for (k = 0; k < nClustersForConsideration; k++)
            {
                l = clusterForConsideration[k];
                qualityValue = edgeWeightPerCluster[l] - network.nodeWeights[j] * clusterWeight[l] * resolution
                               + targetCover.getCoverWeight(j, l) * targetWeight;
                qualityValuePerConsideredCluster[k] = qualityValue;
                if (qualityValue > maxQualityValue)
                {
                    bestCluster = l;
                    maxQualityValue = qualityValue;
                }
                edgeWeightPerCluster[l] = 0;
                addedClusterForConsideration[l] = false;
            }

            // Add the node weight to its new cluster.
            clusterWeight[bestCluster] += network.nodeWeights[j];

            /*
             * Check whether the cluster was empty and if necessary determine
             * the next empty cluster.
             */
            if (nNodesPerCluster[bestCluster] == 0)
            {
                nNodesPerCluster[bestCluster]++;

                // Check to see if there is a smaller empty cluster that
                // we can use instead of the current nextUnusedCluster.
                if (bestCluster <= nextUnusedCluster)
                {
                    nextUnusedCluster = maxNClusters;
                    for (k = bestCluster; k < maxNClusters; k++)
                        if (nNodesPerCluster[k] == 0 && ((k < nodesPerCover.nNodes && nodesPerCover.getNCovers(k) == 0) || k >= nodesPerCover.nNodes))
                        {
                            nextUnusedCluster = k;
                            break;
                        }
                }
            }
            else
            {
                /* Note that we cannot increase the nNodesPerCluster[bestCluster]
                 * outside the if-statement, since it also affects operations
                 * inside the true-block of the if-statement.
                 */
                nNodesPerCluster[bestCluster]++;
            }

            // Mark the currently selected node as stable and remove it from the queue.
            stableNode[j] = true;
            nUnstableNodes--;
            
            /*
             * If the new cluster of the currently selected node is different
             * from the old cluster, some further updating of the clustering
             * statistics is performed. Also, the neighbors of the currently
             * selected node that do not belong to the new cluster are marked as
             * unstable and are added to the queue.
             */
            if (bestCluster != currentCluster)
            {
                improvement += maxQualityValue - currentQualityValue;
                clustering.clusters[j] = bestCluster;
//                System.out.println("Moving node " + j + " to " + bestCluster + " (cover weight " + targetCover.getCoverWeight(j, bestCluster) + "), improves " + (maxQualityValue - currentQualityValue));
                if (bestCluster + 1 > clustering.nClusters)
                    clustering.nClusters = bestCluster + 1;
                for (k = network.firstNeighborIndices[j]; k < network.firstNeighborIndices[j + 1]; k++)
                    if (stableNode[network.neighbors[k]] && (clustering.clusters[network.neighbors[k]] != bestCluster))
                    {
                        stableNode[network.neighbors[k]] = false;
                        nUnstableNodes++;
                        nodeOrder[(i + nUnstableNodes < network.nNodes) ? (i + nUnstableNodes) : (i + nUnstableNodes - network.nNodes)] = network.neighbors[k];
                    }
                update = true;
            }
            i = (i < network.nNodes - 1) ? (i + 1) : 0;
        }
        while (nUnstableNodes > 0);

        return update;
    }
}
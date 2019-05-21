package cwts.networkanalysis;

import java.util.Random;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;

/**
 * Random labeled target cover local merging algorithm.
 *
 * <p>
 * The algorithm always starts from a singleton partition and merges nodes with
 * other cluster. To cluster which to merge with is chosen randomly from all
 * clusters that do not decrease the quality. Clusters leading to larger
 * increases in the quality function have a higher probability to be selected,
 * and the degree to which this is the case is controlled by the randomness. The
 * higher the randomness, the more uniform is this choice. The lower the
 * randomness the more likely the cluster with the highest increase is chosen. A
 * node is merged with another cluster only if both are sufficiently well
 * connected in terms of the quality function to the rest of the network.
 * </p>
 *
 * <p>
 * The local merging is only used for the refinement in the Leiden algorithm.
 * For the refinement, the actual labels of the target cover do not matter. For
 * performance reasons, they are relabeled using {@link
 * Cover#createSubcoverRelabel(int[], int) Cover.createSubcoverRelabel} by the
 * Leiden algorithm, which is not problematic.
 * </p>
 *
 * <p>
 * However, clusters should only be merged whenever they are both sufficiently
 * well connected to the network. For that reason, we must know to what extent a
 * cluster is overlapping with the cover. Since the cover is relabeled, we
 * specify a networkCoverLabel, which indicates the cover of the network as a
 * whole.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 *
 * @see LeidenAlgorithmLabeledTargetCover
 */

public class LocalMergingAlgorithmLabeledTargetCover extends CPMClusteringAlgorithmLabeledTargetCover
{
    /**
     * Default randomness
     */
    public static final double DEFAULT_RANDOMNESS = 1e-2;
    protected double randomness;
    protected Random random;
    protected int networkCoverLabel;

    /**
     * Constructs new labeled target cover local merging algorithm.
     *
     * <p>
     * Uses the default resolution. Uses the default randomness.
     * </p>
     */
    public LocalMergingAlgorithmLabeledTargetCover()
    {
        this(DEFAULT_RESOLUTION, DEFAULT_RANDOMNESS, new Random());
    }

    /**
     * Constructs new labeled target cover local merging algorithm for specified
     * random number generator.
     *
     * <p>
     * Uses the default resolution. Uses the default randomness.
     * </p>
     *
     * @param random Random number generator
     */
    LocalMergingAlgorithmLabeledTargetCover(Random random)
    {
        this(DEFAULT_RESOLUTION, DEFAULT_RANDOMNESS, random);
    }

    /**
     * Constructs new labeled target cover local merging algorithm for specified
     * randomness and random number generator.
     *
     * <p>
     * Uses the default resolution.
     * </p>
     *
     * @param randomness Randomness
     * @param random Random number generator
     */
    public LocalMergingAlgorithmLabeledTargetCover(double randomness, Random random)
    {
        this(DEFAULT_RESOLUTION, randomness, random);
    }

    /**
     * Constructs new labeled target cover local merging algorithm for the
     * specified resolution, randomness and random number generator.
     *
     * @param resolution Resolution parameter
     * @param randomness Randomness
     * @param random Random number generator
     */
    public LocalMergingAlgorithmLabeledTargetCover(double resolution, double randomness, Random random)
    {
        super(resolution);
        this.randomness = randomness;
        this.random = random;
    }

    /**
     * Gets the network cover label.
     *
     * @return Network cover label
     */
    public double getNetworkCoverLabel()
    {
        return this.networkCoverLabel;
    }

    /**
     * Gets randomness
     *
     * @return Randomness
     */
    public double getRandomness()
    {
        return randomness;
    }

    /**
     * Gets random number generator
     *
     * @return Random number generator
     */
    public Random getRandom()
    {
        return random;
    }


    /**
     * Sets the network cover label.
     *
     * @param networkCoverLabel Network cover label
     */
    public void setNetworkCoverLabel(int networkCoverLabel)
    {
        this.networkCoverLabel = networkCoverLabel;
    }

    /**
     * Sets randomness
     *
     * @param randomness Randomness
     */
    public void setRandomness(double randomness)
    {
        this.randomness = randomness;
    }

    /**
     * Sets random number generator
     *
     * @param random Random number generator
     *
     */
    public void setRandom(Random random)
    {
        this.random = random;
    }

    /**
     * Runs local merging algorithm.
     *
     * <p>
     * The algorithm always starts from a singleton partition and merges nodes
     * with other cluster. To cluster which to merge with is chosen randomly
     * from all clusters that do not decrease the quality. Clusters leading to
     * larger increases in the quality function have a higher probability to be
     * selected, and the degree to which this is the case is controlled by the
     * randomness. The higher the randomness, the more uniform is this choice.
     * The lower the randomness the more likely the cluster with the highest
     * increase is chosen. A node is merged with another cluster only if both
     * are sufficiently well connected in terms of the quality function to the
     * rest of the network.
     * </p>
     *
     * <p>
     * The local merging is only used for the refinement phase in the Leiden
     * algorithm. For the refinement, the actual labels of the target cover do
     * not matter. For performance reasons, they are relabeled using {@link
     * Cover#createSubcoverRelabel(int[], int) Cover.createSubcoverRelabel} by
     * the Leiden algorithm, which is not problematic.
     * </p>
     *
     * <p>
     * However, clusters should only be merged whenever they are both
     * sufficiently well connected to the network. For that reason, we must know
     * to what extent a cluster is overlapping with the cover. Since the cover
     * is relabeled, we specify a networkCoverLabel, which indicates the cover
     * of the network as a whole.
     * </p>
     *
     * <p>
     * Empty clusters are always removed. Since the actual labels do not
     * matter in the refinement phase, we can safely remove the empty clusters.
     * </p>
     *
     * @return Boolean indicating whether the clustering has changed
     */
    public Clustering findClustering(Network network)
    {
        boolean update;
        boolean[] nonSingletonCluster, mayMerge, addedClusterForConsideration;
        double maxQualityValueIncrement, qualityValueIncrement, r, totalNodeWeight, totalTransformedQualityValueIncrement;
        double[] clusterCoverWeightNetworkCoverLabel, clusterWeight, cumTransformedQualityValueIncrementPerCluster, edgeWeightPerCluster, externalEdgeWeightPerCluster;
        int bestCluster, chosenCluster, currentCluster, i, j, k, l, max_idx, mid_idx, min_idx, nNeighboringClusters;
        int[] neighboringCluster, nodeOrder;

        // Start from singleton clustering.
        Clustering clustering = new Clustering(network.nNodes);

        if (network.nNodes == 1)
            return clustering;

        update = false;

        int maxNClusters = Math.max(network.nNodes, targetCover.nCovers) + network.nNodes;
        clustering.nClusters = maxNClusters;
        
        totalNodeWeight = network.getTotalNodeWeight();
        clusterWeight = new double[maxNClusters];
        nonSingletonCluster = new boolean[maxNClusters];
        externalEdgeWeightPerCluster = new double[maxNClusters];

        nodeOrder = cwts.util.Arrays.generateRandomPermutation(network.nNodes, random);

        clusterCoverWeightNetworkCoverLabel = new double[maxNClusters];
        edgeWeightPerCluster = new double[maxNClusters];
        neighboringCluster = new int[maxNClusters];
        cumTransformedQualityValueIncrementPerCluster = new double[maxNClusters];
        addedClusterForConsideration = new boolean[maxNClusters];

        for (i = 0; i < network.nNodes; i++)
        {
            clusterWeight[i] = network.nodeWeights[i];
            externalEdgeWeightPerCluster[i] = network.getTotalEdgeWeight(i);
            clusterCoverWeightNetworkCoverLabel[i] = targetCover.getCoverWeight(i, networkCoverLabel);
        }

        // Determine whether cluster may merge
        mayMerge = new boolean[maxNClusters]; // Initially set to false
        for (i = 0; i < network.nNodes; i++)
            if (externalEdgeWeightPerCluster[i] >= clusterWeight[i] * (totalNodeWeight - clusterWeight[i]) * resolution
                                                   - clusterCoverWeightNetworkCoverLabel[i] * targetWeight)
                mayMerge[i] = true;
        for (; i < maxNClusters; i++)
                mayMerge[i] = true;
        
        for (i = 0; i < network.nNodes; i++)
        {
            j = nodeOrder[i];
            currentCluster = clustering.clusters[j];

            /*
             * Only nodes belonging to singleton clusters can be moved to a
             * different cluster. This guarantees that clusters will never be
             * split up. Additionally, only nodes that are well connected with
             * the rest of the network are considered for moving.
             */
            if (!nonSingletonCluster[currentCluster] && mayMerge[currentCluster])
            {
                /*
                 * Only nodes belonging to singleton clusters can be moved to a
                 * different cluster. This guarantees that clusters will never
                 * be split up. Additionally, only nodes that are well connected
                 * with the rest of the network are considered for moving.
                 */
                clusterWeight[currentCluster] = 0;
                externalEdgeWeightPerCluster[currentCluster] = 0;
                clusterCoverWeightNetworkCoverLabel[currentCluster] = 0;

                // Identify the neighboring clusters of the currently selected node, that is, the
                // clusters with which the currently selected node is connected. The old cluster of
                // the currently selected node is also included in the set of neighboring clusters.
                // In this way, it is always possible that the currently selected node will be
                // moved back to its old cluster.
                neighboringCluster[0] = currentCluster;
                addedClusterForConsideration[neighboringCluster[0]] = true;
                nNeighboringClusters = 1;
                for (k = network.firstNeighborIndices[j]; k < network.firstNeighborIndices[j + 1]; k++)
                {
                    l = clustering.clusters[network.neighbors[k]];
                    if (!addedClusterForConsideration[l])
                    {
                        neighboringCluster[nNeighboringClusters] = l;
                        addedClusterForConsideration[l] = true;
                        nNeighboringClusters++;
                    }
                    edgeWeightPerCluster[l] += network.edgeWeights[k];
                }

                // Also add covers of node
                for(Int2DoubleMap.Entry entry : targetCover.getCoverIterator(j))
                {
                    l = entry.getIntKey();
                    if (!addedClusterForConsideration[l])
                    {
                        neighboringCluster[nNeighboringClusters] = l;
                        addedClusterForConsideration[l] = true;
                        nNeighboringClusters++;
                    }
                }

                /*
                 * For each cluster for consideration, determine whether the
                 * neighboring cluster is well connected with the rest of the
                 * network. For each cluster for consideration that is well
                 * connected, calculate the increment of the quality function
                 * obtained by moving the currently selected node to the
                 * cluster. For each cluster for consideration for which the
                 * increment is non-negative, calculate a transformed increment
                 * that will determine the probability with which the currently
                 * selected node is moved to the neighboring cluster.
                 */
                bestCluster = currentCluster;
                maxQualityValueIncrement = 0;
                totalTransformedQualityValueIncrement = 0;
                for (k = 0; k < nNeighboringClusters; k++)
                {
                    l = neighboringCluster[k];

                    if (mayMerge[l])
                    {
                        qualityValueIncrement = edgeWeightPerCluster[l] - network.nodeWeights[j] * clusterWeight[l] * resolution
                                                + targetCover.getCoverWeight(j, l) * targetWeight;

                        if (qualityValueIncrement > maxQualityValueIncrement)
                        {
                            bestCluster = l;
                            maxQualityValueIncrement = qualityValueIncrement;
                        }

                        if (qualityValueIncrement >= 0)
                            totalTransformedQualityValueIncrement += cwts.util.FastMath.fastExp(qualityValueIncrement / randomness);
                    }

                    cumTransformedQualityValueIncrementPerCluster[k] = totalTransformedQualityValueIncrement;

                    edgeWeightPerCluster[l] = 0;
                    addedClusterForConsideration[l] = false;
                }

                /*
                 * Determine the neighboring cluster to which the currently
                 * selected node will be moved.
                 */
                if (maxQualityValueIncrement > 0)
                {
                    if (totalTransformedQualityValueIncrement < Double.POSITIVE_INFINITY)
                    {
                        r = totalTransformedQualityValueIncrement * random.nextDouble();
                        min_idx = -1;
                        max_idx = nNeighboringClusters + 1;
                        while (min_idx < max_idx - 1)
                        {
                            mid_idx = (min_idx + max_idx) / 2;
                            if (cumTransformedQualityValueIncrementPerCluster[mid_idx] >= r)
                                max_idx = mid_idx;
                            else
                                min_idx = mid_idx;
                        }
                        chosenCluster = neighboringCluster[max_idx];
                    }
                    else
                        chosenCluster = bestCluster;
                }
                else
                    chosenCluster = currentCluster;

                // Add the node weight to its new cluster.
                clusterWeight[chosenCluster] += network.nodeWeights[j];
                clusterCoverWeightNetworkCoverLabel[chosenCluster] += targetCover.getCoverWeight(j, networkCoverLabel);

                /*
                 * Recalculate the number of edges of each cluster to the rest
                 * of the network.
                 */
                for (k = network.firstNeighborIndices[j]; k < network.firstNeighborIndices[j + 1]; k++)
                    if (clustering.clusters[network.neighbors[k]] == chosenCluster)
                        externalEdgeWeightPerCluster[chosenCluster] -= network.edgeWeights[k];
                    else
                        externalEdgeWeightPerCluster[chosenCluster] += network.edgeWeights[k];

                /*
                 * Recalculate wether the chosen cluster may still merge, and
                 * indicate that the chosen cluster is now not a singleton
                 * anymore.
                 */
                if (chosenCluster != j)
                {
                    clustering.clusters[j] = chosenCluster;
                    if (externalEdgeWeightPerCluster[chosenCluster] < clusterWeight[chosenCluster] * (totalNodeWeight - clusterWeight[chosenCluster]) * resolution
                                                                      - clusterCoverWeightNetworkCoverLabel[chosenCluster] * targetWeight)
                        mayMerge[chosenCluster] = false;

                    nonSingletonCluster[chosenCluster] = true;
                    update = true;
                }
            }
        }

        clustering.removeEmptyClusters();

        return clustering;
    }
}
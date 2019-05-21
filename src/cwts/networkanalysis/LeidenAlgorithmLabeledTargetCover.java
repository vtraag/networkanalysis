package cwts.networkanalysis;

import java.util.Arrays;
import java.util.Random;

/**
 * Labeled target cover Leiden algorithm.
 *
 * <p>
 * The Leiden algorithm consists of three phases: (1) local moving of nodes, (2)
 * refinement of the partition and (3) aggregation of the network based on the
 * refined partition, using the non-refined partition to create an initial
 * partition for the aggregate network. These phases are repeated until nothing
 * can be aggregated further. The local moving is done by the
 * FastLocalMovingAlgorithm. The refinement is done by the
 * RandomLocalMergingAlgorithm. All these algorithms are stochastic. If no
 * explicit random object is provided, they will be initialized using a random
 * seed.
 * </p>
 *
 * <p>
 * The resolution parameter of any local moving algorithm is always set to the
 * resolution parameter of the Leiden algorithm.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 *
 * @see FastLocalMovingAlgorithmLabeledTargetCover
 * @see RandomLocalMergingAlgorithmLabeledTargetCover
 */
public class LeidenAlgorithmLabeledTargetCover extends IterativeCPMClusteringAlgorithmLabeledTargetCover
{
    /**
     * Default randomness parameter.
     */
    public static final double DEFAULT_RANDOMNESS = LocalMergingAlgorithm.DEFAULT_RANDOMNESS;

    /**
     * Randomness parameter.
     */
    protected double randomness;

    /**
     * Random number generator.
     */
    protected Random random;
    
    /**
     * Local moving algorithm.
     */
    protected IncrementalCPMClusteringAlgorithmLabeledTargetCover localMovingAlgorithm;
    
    /**
     * Constructs new labeled target cover Leiden algorithm for specified
     * network, resolution parameter, labeled target cover, target weight and
     * random number generator.
     *
     * <p>
     * Uses a singleton clustering. Uses the default number of iterations. Uses
     * the fast target cover local moving algorithm overall. Uses the random
     * target cover local merging algorithm per cluster. Both local moving
     * algorithms are initialized with the indicated random number generator.
     * </p>
     *
     * @param resolution Resolution parameter
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     * @param randomness  Randomness parameter
     * @param random Random number generator
     */
    public LeidenAlgorithmLabeledTargetCover(double resolution, Cover targetCover, double targetWeight)
    {
        this(resolution, DEFAULT_N_ITERATIONS, new FastLocalMovingAlgorithmLabeledTargetCover(new Random()), targetCover, targetWeight, DEFAULT_RANDOMNESS, new Random());
    }
    
    /**
     * Constructs new labeled target cover Leiden algorithm for specified
     * network, resolution parameter, labeled target cover, target weight and
     * random number generator.
     *
     * <p>
     * Uses a singleton clustering. Uses the default number of iterations. Uses
     * the fast target cover local moving algorithm overall. Uses the random
     * target cover local merging algorithm per cluster. Both local moving
     * algorithms are initialized with the indicated random number generator.
     * </p>
     *
     * @param resolution Resolution parameter
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     * @param randomness  Randomness parameter
     * @param random Random number generator
     */
    public LeidenAlgorithmLabeledTargetCover(double resolution, Cover targetCover, double targetWeight, double randomness, Random random)
    {
        this(resolution, DEFAULT_N_ITERATIONS, new FastLocalMovingAlgorithmLabeledTargetCover(random), targetCover, targetWeight, randomness, random);
    }

    /**
     * Constructs new labeled target cover Leiden algorithm for specified
     * network, clustering, resolution parameter, number of iterations, overall
     * local moving algorithm, local moving algorithm per cluster, labeled
     * target cover and target weight.
     *
     * @param network Network
     * @param clustering Clustering
     * @param resolution Resolution parameter
     * @param nIterations Number of iterations
     * @param targetCover Target cover
     * @param targetWeight Target weight
     */
    public LeidenAlgorithmLabeledTargetCover(double resolution, int nIterations, Cover targetCover, double targetWeight, double randomness, Random random)
    {
        this(resolution, nIterations, new FastLocalMovingAlgorithmLabeledTargetCover(random), targetCover, targetWeight, randomness, random);
    }
    
    /**
     * Constructs new labeled target cover Leiden algorithm for specified
     * network, clustering, resolution parameter, number of iterations, overall
     * local moving algorithm, local moving algorithm per cluster, labeled
     * target cover and target weight.
     *
     * @param network Network
     * @param clustering Clustering
     * @param resolution Resolution parameter
     * @param nIterations Number of iterations
     * @param localMovingAlgorithm Labeled target cover local moving algorithm
     * @param targetCover Target cover
     * @param targetWeight Target weight
     */
    public LeidenAlgorithmLabeledTargetCover(double resolution, int nIterations, IncrementalCPMClusteringAlgorithmLabeledTargetCover localMovingAlgorithm, Cover targetCover, double targetWeight, double randomness, Random random)
    {
        super(resolution, nIterations, targetCover, targetWeight);
        
        this.randomness = randomness;
        this.random = random;
        setLocalMovingAlgorithm(localMovingAlgorithm);
    }

    /**
     * Gets the labeled target cover local moving algorithm.
     *
     * @return Overall local moving algorithm
     */
    public IncrementalCPMClusteringAlgorithmLabeledTargetCover getLocalMovingAlgorithm()
    {
        return localMovingAlgorithm;
    }

    /**
     * Sets the overall local labeled target cover moving algorithm.
     *
     * @param localMovingAlgorithm Overall local moving algorithm
     */
    public void setLocalMovingAlgorithm(IncrementalCPMClusteringAlgorithmLabeledTargetCover localMovingAlgorithm)
    {
        this.localMovingAlgorithm = localMovingAlgorithm;
        this.localMovingAlgorithm.resolution = resolution;
        this.localMovingAlgorithm.targetCover = targetCover;
        this.localMovingAlgorithm.targetWeight = targetWeight;
    }

    /**
     * Runs single iteration of the Leiden algorithm.
     *
     * <p>
     * The Leiden algorithm consists of three phases: (1) local moving of nodes,
     * (2) refinement of the partition and (3) aggregation of the network based
     * on the refined partition, using the non-refined partition to create an
     * initial partition for the aggregate network. These phases are repeated
     * until nothing can be aggregated further.
     * </p>
     *
     * @return Boolean indicating whether the clustering has changed
     */
    @Override
    protected boolean improveClusteringOneIteration(Network network, Clustering clustering)
    {
        boolean update;
        double[] subnetworkEdgeWeight;
        int i, j, k;
        int[] nNodesPerClusterReducedNetwork, subnetworkNeighbor, subnetworkNode;
        int[][] nodePerCluster;
        LeidenAlgorithmLabeledTargetCover leidenAlgorithm;
        Cover subcover;
        Network subnetwork;
        Clustering oldClustering;
        
        // Update the partition by moving individual nodes between clusters.
        localMovingAlgorithm.targetCover = targetCover;
        localMovingAlgorithm.targetWeight = targetWeight;
        update = localMovingAlgorithm.improveClustering(network, clustering);

        /*
         * Relabel clusters without overlap, and remove empty clusters that are
         * larger than the number of covers.
         */
        clustering.relabelClustersWithoutOverlap(targetCover, targetCover.getNCovers());
        clustering.removeEmptyClustersLargerThan(targetCover.getNCovers());

        /*
         * If all clusters are nodes, and nothing is aggregated in the local
         * move, we can stop. We have to use the number of non-empty clusters,
         * otherwise we might not converge, because the labeled target cover may
         * need more clusters than nodes.
         */
        if (clustering.getNNonEmptyClusters() < network.nNodes) 
        {
            /*
             * In order to save memory, we only generate each subcluster while
             * iterating. The following variables are necessary to do so.
             */
            oldClustering = new Clustering(clustering.clusters);
            nodePerCluster = oldClustering.getNodesPerCluster();
            subnetworkNode = new int[network.nNodes];
            subnetworkNeighbor = new int[network.nEdges];
            subnetworkEdgeWeight = new double[network.nEdges];
            
            int nSubnetworks = oldClustering.nClusters;

            LocalMergingAlgorithmLabeledTargetCover localMergingAlgorithm = new LocalMergingAlgorithmLabeledTargetCover(resolution, randomness, random);

            /*
             * Refine the partition by iterating over the clusters and by trying
             * to split up each cluster into multiple clusters.
             */
            clustering.nClusters = 0;
            nNodesPerClusterReducedNetwork = new int[nSubnetworks];
            for (i = 0; i < nSubnetworks; i++)
            {
                subnetwork = network.createSubnetwork(oldClustering, i, nodePerCluster[i], subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight);

                /*
                 * Relabel the subcover, since the actual labels are irrelevant
                 * for the refinement. They are only used for aggregation, not
                 * for the actual assignment to labeled clusters.
                 */
                subcover = targetCover.createSubcoverRelabel(nodePerCluster[i], i);

                // Run refinement on subnetwork cluster
                localMergingAlgorithm.targetCover = subcover;
                localMergingAlgorithm.setNetworkCoverLabel(0);
                Clustering clusteringSubnetwork = localMergingAlgorithm.findClustering(subnetwork);

                /*
                 * Update the overall clustering that will be used for
                 * refinement.
                 */
                for (j = 0; j < subnetwork.nNodes; j++)
                    clustering.clusters[nodePerCluster[i][j]] = clustering.nClusters + clusteringSubnetwork.clusters[j];
                clustering.nClusters += clusteringSubnetwork.nClusters;
                nNodesPerClusterReducedNetwork[i] = clusteringSubnetwork.nClusters;
            }

            /*
             * If nothing was aggregated, aggregate complete clusters. It is
             * possible that nothing, or only very little was aggregated in the
             * local moving, especially because of the target cover. We
             * aggregate the complete clusters in that case to make sure that
             * the aggregation is sufficiently fast.
             */
            if (clustering.nClusters == network.nNodes)
            {
                int c = 0;
                for (i = 0; i < nSubnetworks; i++)
                {
                    for (j = 0; j < nodePerCluster[i].length; j++)    
                        clustering.clusters[nodePerCluster[i][j]] = c;
                    if (nodePerCluster[i].length > 0)
                    {
                        c += 1;
                        nNodesPerClusterReducedNetwork[i] = 1; 
                    }
                    else
                        nNodesPerClusterReducedNetwork[i] = 0; 
                }       
                clustering.nClusters = c;
            }

            /*
             * Create new Leiden algorithm on aggregate network based on the
             * refined partition of the non-aggregate network.
             */
            Network reducedNetwork = network.createReducedNetwork(clustering);
            Cover reducedCover = targetCover.createReducedCover(clustering);            
            leidenAlgorithm = new LeidenAlgorithmLabeledTargetCover(resolution, 1,
                                                             localMovingAlgorithm, 
                                                             reducedCover, targetWeight,
                                                             randomness, random);
            /*
             * Create an initial clustering for the aggregate network based on
             * the non-refined clustering of the non-aggregate network.
             */
            int[] clustersReducedNetwork = new int[clustering.nClusters];
            i = 0;
            for (j = 0; j < nNodesPerClusterReducedNetwork.length; j++)
            {
                Arrays.fill(clustersReducedNetwork, i, i + nNodesPerClusterReducedNetwork[j], j);
                i += nNodesPerClusterReducedNetwork[j];
            }
            Clustering clusteringReducedNetwork = new Clustering(clustersReducedNetwork);
            
            /*
             * Recursively run the algorithm to the aggregate network. The
             * algorithm starts from the initial partition adapted above.
             */
            update |= leidenAlgorithm.improveClustering(reducedNetwork, clusteringReducedNetwork);

            /*
             * Update the partition of the non-aggregate network so that it
             * coincides with the final partition obtained for the aggregate
             * network.
             */
            clustering.mergeClusters(clusteringReducedNetwork);
        }

        /*
         * Relabel clusters without overlap, and remove empty clusters that are
         * larger than the number of covers.
         */
        clustering.relabelClustersWithoutOverlap(targetCover, targetCover.getNCovers());
        clustering.removeEmptyClustersLargerThan(targetCover.getNCovers());        
        return update;
    }
}
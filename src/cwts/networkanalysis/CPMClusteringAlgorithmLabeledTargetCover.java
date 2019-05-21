package cwts.networkanalysis;

/**
 * Base class for CPM based algorithms, taking into consideration a
 * labeled target cover.
 *
 * <p>
 * Provides the implementation of the quality function, and enables the
 * specification of a resolution parameter. If no resolution is specified, the
 * default resolution is 1.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 *
 */
public abstract class CPMClusteringAlgorithmLabeledTargetCover extends ClusteringAlgorithmTargetCover
{

    /**
     * Default resolution parameter.
     */
    public static final double DEFAULT_RESOLUTION = 1;
    protected double resolution;

    /**
     * Constructs new labeled target cover CPM clustering algorithm.
     *
     * <p>
     * Uses the default resolution parameter.
     * </p>
     */
    public CPMClusteringAlgorithmLabeledTargetCover()
    {
        this(DEFAULT_RESOLUTION);
    }

    /**
     * Constructs new labeled target cover CPM clustering algorithm for
     * specified resolution parameter.
     *
     * @param resolution Resolution parameter
     */
    public CPMClusteringAlgorithmLabeledTargetCover(double resolution)
    {
        super();
        this.resolution = resolution;
    }

    /**
     * Constructs new labeled target cover CPM clustering algorithm for
     * specified network, resolution, labeled target cover and target weight.
     *
     * <p>
     * Uses singleton clustering.
     * </p>
     * @param resolution Resolution parameter
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public CPMClusteringAlgorithmLabeledTargetCover(double resolution,  Cover targetCover, double targetWeight)
    {
        super(targetCover, targetWeight);
        this.resolution = resolution;
    }

    /**
     * Gets resolution parameter.
     *
     * @return Resolution parameter
     */
    public double getResolution()
    {
        return resolution;
    }

    /**
     * Sets resolution parameter.
     *
     * @param resolution Resolution parameter
     */
    public void setResolution(double resolution)
    {
        this.resolution = resolution;
    }

    /**
     * Calculates quality of a partition according to CPM using a labeled
     * target cover.
     *
     * <p>
     * Labeled target cover CPM is defined as follows:
     * </p>
     *
     * <blockquote>
     * 1/(2m) (sum_{ij} [A_{ij} - gamma n_i n_j]d(s_i, s_j) + lambda sum_i S_{i s_i}).
     * </blockquote>
     *
     * <p>
     * Here m is the sum of all link weights in the network, A_{ij} is the
     * weight of the link between i and j, which equals zero if there is no
     * link, gamma is the resolution parameter, n_i is the node size of node i,
     * d(s_i, s_j) = 1 if s_i = s_j and equals zero otherwise, where s_i is the
     * cluster of node i, S_{id} denotes the weight of cover d for node i in the
     * labeled target cover, so that S_{i s_i} denotes the weight of cover s_i for node
     * i, and lambda denotes the target weight.
     * </p>
     *
     * @return Labeled target cover CPM
     */
    public double calcQuality(Network network, Clustering clustering)
    {
        double qualityValue;
        double[] clusterWeight;
        int i, j, k;
        qualityValue = 0;

        // Add all internal edge weights to the quality value
        for (i = 0; i < network.nNodes; i++)
        {
            j = clustering.clusters[i];
            for (k = network.firstNeighborIndices[i]; k < network.firstNeighborIndices[i + 1]; k++)
                if (clustering.clusters[network.neighbors[k]] == j)
                    qualityValue += network.edgeWeights[k];
        }

        // Add self loops to the quality value
        qualityValue += network.totalEdgeWeightSelfLinks;

        /*
         * Note that sum_ij n_i n_j d(s_i, s_j) can be written as sum_c n_c^2
         * where n_c = sum_i n_i d(s_i, c). We here calculate the squares of
         * cluster sizes, multiply by the resolution, and subtract that from the
         * quality value.
         */
        clusterWeight = new double[clustering.nClusters];
        double clusterWeightCosts = 0.0;
        for (i = 0; i < network.nNodes; i++)
            clusterWeight[clustering.clusters[i]] += network.nodeWeights[i];
        for (i = 0; i < clustering.nClusters; i++)
            clusterWeightCosts += clusterWeight[i] * clusterWeight[i] * resolution;
        qualityValue -= clusterWeightCosts;

        /*
         * Note that sum_i S_{i s_i} can be written as
         * sum_c (sum_i S_{i s_i}\delta(s_i, c)).
         * We here calculate cover weights for all clusters, and sum the cover
         * weight of cluster i with cover i, multiply by the target weight,
         * and add that to the quality value.
         */
        double coverBenefits = 0.0;
        Cover reducedCover = this.targetCover.createReducedCover(clustering);
        for (i = 0; i < reducedCover.nNodes; i++)
        {
            // Only count to what extent the actual assignment overlaps with their labeled target cover
            coverBenefits += reducedCover.getCoverWeight(i, i) * this.targetWeight;
        }
        qualityValue += coverBenefits;

        // Divide by the total edge weight to normalize
        qualityValue /= 2 * network.getTotalEdgeWeight() + network.totalEdgeWeightSelfLinks;
        return qualityValue;
    }
}
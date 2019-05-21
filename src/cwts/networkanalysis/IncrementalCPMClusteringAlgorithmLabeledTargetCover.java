package cwts.networkanalysis;

/**
 * Abstract base class for incremental clustering algorithms that use the CPM
 * quality function.
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 */
public abstract class IncrementalCPMClusteringAlgorithmLabeledTargetCover extends CPMClusteringAlgorithmLabeledTargetCover implements IncrementalClusteringAlgorithm
{

    /**
     * Constructs new labeled target cover CPM clustering algorithm for
     * specified resolution.
     *
     * @param resolution Resolution parameter
     */
    public IncrementalCPMClusteringAlgorithmLabeledTargetCover(double resolution)
    {
        super(resolution);
    }

    /**
     * Constructs new labeled target cover CPM clustering algorithm for
     * specified resolution, labeled target cover and target weight.
     *
     * @param resolution Resolution parameter
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public IncrementalCPMClusteringAlgorithmLabeledTargetCover(double resolution,  Cover targetCover, double targetWeight)
    {
        super(resolution, targetCover, targetWeight);
    }

    /**
     * Finds a clustering of the nodes in a network.
     *
     * <p>
     * The clustering is obtained by calling {@link #improveClustering(Network
     * network, Clustering clustering)} and by providing a singleton clustering
     * as input to this method.
     * </p>
     *
     * @param network Network
     *
     * @return Clustering
     */
    public Clustering findClustering(Network network)
    {
        Clustering clustering;

        clustering = new Clustering(network.getNNodes());
        improveClustering(network, clustering);
        return clustering;
    }
}

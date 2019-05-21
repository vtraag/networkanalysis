package cwts.networkanalysis;
/**
 * Base class for modularity based clustering algorithm, taking into
 * consideration a labeled target cover, that is iterated multiple times.
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 */
public abstract class IterativeCPMClusteringAlgorithmLabeledTargetCover extends IncrementalCPMClusteringAlgorithmLabeledTargetCover
{

    /**
     * Default number of iterations (i.e. iterate as long as there are changes).
     */
    public static final int DEFAULT_N_ITERATIONS = 1;
    protected int nIterations;

    /**
     * Constructs new iterated labeled target cover modularity algorithm for
     * specified resolution parameter, target cover and target weight using
     * the default number of iterations.
     *
     * @param resolution Resolution parameter
     */
    public IterativeCPMClusteringAlgorithmLabeledTargetCover(double resolution)
    {
        super(resolution);
    }

    /**
     * Constructs new iterated labeled target cover modularity algorithm for
     * specified resolution parameter, target cover and target weight using
     * the default number of iterations.
     *
     * @param resolution Resolution parameter
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public IterativeCPMClusteringAlgorithmLabeledTargetCover(double resolution, Cover targetCover, double targetWeight)
    {
        this(resolution, DEFAULT_N_ITERATIONS, targetCover, targetWeight);
    }
    
    /**
     * Constructs new iterated labeled target cover modularity algorithm for
     * specified resolution parameter, number of iterations,
     * target cover and target weight.
     *
     * @param resolution Resolution parameter
     * @param nIterations Number of iterations
     * @param targetCover Labeled target cover
     * @param targetWeight Target weight
     */
    public IterativeCPMClusteringAlgorithmLabeledTargetCover(double resolution, int nIterations, Cover targetCover, double targetWeight)
    {
        super(resolution, targetCover, targetWeight);
        this.nIterations = nIterations;
    }

    /**
     * Gets number of iterations.
     *
     * @return Number of iterations
     */
    public int getNIterations()
    {
        return nIterations;
    }

    /**
     * Sets the number of iterations.
     *
     * @param nIterations Number of iterations
     */
    public void setNIterations(int nIterations)
    {
        this.nIterations = nIterations;
    }

    /**
     * Improves a clustering of the nodes in a network.
     *
     * <p>
     * If the number of iterations {@code nIterations} is positive, the
     * clustering is improved by making {@code nIterations} calls to {@link
     * #improveClusteringOneIteration(Network network, Clustering clustering)}.
     * If {@code nIterations} equals 0, calls to {@link
     * #improveClusteringOneIteration(Network network, Clustering clustering)}
     * continue to be made until there has been a call that did not result in
     * an improvement of the clustering.
     * </p>
     *
     * @param network    Network
     * @param clustering Clustering
     *
     * @return Boolean indicating whether the clustering has been improved
     */
    public boolean improveClustering(Network network, Clustering clustering)
    {
        boolean update;
        int i;

        update = false;
        if (nIterations > 0)
            for (i = 0; i < nIterations; i++)
                update |= improveClusteringOneIteration(network, clustering);
        else
            while (improveClusteringOneIteration(network, clustering))
                update = true;
        return update;
    }

    /**
     * Improves a clustering by performing one iteration of an iterative
     * clustering algorithm.
     *
     * @param network    Network
     * @param clustering Clustering
     *
     * @return Boolean indicating whether the clustering has been improved
     */
    protected abstract boolean improveClusteringOneIteration(Network network, Clustering clustering);
}
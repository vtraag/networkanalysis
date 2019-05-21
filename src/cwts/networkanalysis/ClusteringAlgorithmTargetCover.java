package cwts.networkanalysis;

/**
 * Base class for clustering algorithm, taking into consideration a target cover.
 *
 * <p>
 * A target cover is understood to refer to a (possibly overlapping) clustering
 * that is taken into consideration when clustering a network. It is assumed
 * that if the target weight is very high, then the clustering returned is
 * highly similar to the target cover.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag */
public abstract class ClusteringAlgorithmTargetCover implements ClusteringAlgorithm
{

    protected double targetWeight;
    protected Cover targetCover;

    /**
     * Default target weight.
     */
    public static final double DEFAULT_TARGETWEIGHT = 0;

    /**
     * Constructs target cover clustering algorithm without any network,
     * clustering or target cover.
     */
    public ClusteringAlgorithmTargetCover() {}

    /**
     * Constructs target cover clustering algorithm for specified network,
     * clustering, target cover and target weight.
     
     * @param targetCover Target cover
     * @param targetWeight Target weight
     */
    public ClusteringAlgorithmTargetCover(Cover targetCover, double targetWeight)
    {
        this.targetCover = targetCover;
        this.targetWeight = targetWeight;
    }

    /**
     * Gets target cover.
     *
     * @return Target cover
     */
    public Cover getTargetCover()
    {
        return targetCover;
    }

    /**
     * Sets target cover.
     *
     * @param targetCover Target cover
     */
    public void setTargetCover(Cover targetCover)
    {
        this.targetCover = targetCover;
    }

    /**
     * Gets target weight.
     *
     * @return Target weight
     */
    public double getTargetWeight()
    {
        return targetWeight;
    }

    /**
     * Sets target weight.
     *
     * @param targetWeight Target weight
     */
    public void setTargetWeight(double targetWeight)
    {
        this.targetWeight = targetWeight;
    }
}
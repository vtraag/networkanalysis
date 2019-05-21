package cwts.networkanalysis;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;

/**
 * Cover of nodes of a graph.
 *
 * <p>
 * Each node {@code 0, 1, ...., nNodes} is assigned to zero or more covers
 * {@code 0, 1, ..., nCovers} with a certain weight.
 * </p>
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Vincent Traag
 */
public class Cover implements Cloneable, Serializable
{

    private static final long serialVersionUID = 1;
    int nNodes;
    int nCovers;
    Int2DoubleOpenHashMap[] coverMap;

    /**
     * Loads Cover object from file.
     *
     * <p>
     * Should be saved using the save function.
     * </p>
     *
     * @param fileName File in which Cover object was stored
     * @return Loaded Clustering
     * @throws ClassNotFoundException Class Cover not found
     * @throws IOException File not found
     * @see #save(String filename)
     */
    public static Cover load(String fileName) throws ClassNotFoundException, IOException
    {
        Cover covering;
        ObjectInputStream objectInputStream;
        objectInputStream = new ObjectInputStream(new FileInputStream(fileName));
        covering = (Cover) objectInputStream.readObject();
        objectInputStream.close();
        return covering;
    }

    /**
     * Constructs singleton cover for specified number of nodes.
     *
     * @param nNodes Number of nodes
     */
    public Cover(int nNodes)
    {
        this(nNodes, true);
    }

    /**
     * Constructs cover for specified number of nodes.
     *
     * <p>
     * If {@code createSingletonCover = true} then a singleton cover will be
     * created.
     * </p>
     *
     * @param nNodes Number of nodes
     * @param createSingletonCover Create singleton cover
     */
    protected Cover(int nNodes, boolean createSingletonCover)
    {
        this.nNodes = nNodes;
        if (createSingletonCover)
        {
            this.coverMap = new Int2DoubleOpenHashMap[nNodes];
            for (int i = 0; i < nNodes; i++)
            {
                coverMap[i] = new Int2DoubleOpenHashMap();
            }
        }
    }

    /**
     * Constructs cover for specified cover map.
     *
     * <p>
     * The cover map is a an array of int-double hashmaps. That is,
     * {@code coverMap[i].get(c)} is the weight with which node i belongs to
     * cover c.
     * </p>
     *
     * @param coverMap Cover map
     */
    public Cover(Int2DoubleOpenHashMap[] coverMap)
    {
        this.nNodes = coverMap.length;
        this.coverMap = coverMap.clone();
        this.recalcNCovers();
    }

    /**
     * Construct cover from specified cover memberships and cover weights.
     *
     * <p>
     * {@code cover[i].length} is the number of covers in which node i is
     * contained, which should equal {@code coverWeight[i].length}. Cover
     * {@code cover[i][k]} is the k-th cover and node i belongs to the k-th
     * cover with weight {@code coverWeight[i][k]}.
     * </p>
     *
     * @param cover Cover membership
     * @param coverWeight Cover weights
     */
    public Cover(int[][] cover, double[][] coverWeight)
    {
        nNodes = cover.length;
        nCovers = 0;
        this.coverMap = new Int2DoubleOpenHashMap[nNodes];
        for (int i = 0; i < cover.length; i++)
        {
            coverMap[i] = new Int2DoubleOpenHashMap(cover[i].length);
            for (int c = 0; c < cover[i].length; c++)
            {
                coverMap[i].put(cover[i][c], coverWeight[i][c]);
                if (cover[i][c] + 1 > nCovers)
                {
                    nCovers = cover[i][c] + 1;
                }
            }
            coverMap[i].trim();
        }
    }

    /**
     * Construct cover from specified clustering.
     *
     * <p>
     * Each node i is assigned to a single cover, as indicated by its cluster.
     * If the cluster of a node is negative, the node is not assigned to any
     * cover.
     * </p>
     *
     * @param clustering Clustering
     */
    public Cover(Clustering clustering)
    {
        this(clustering, clustering.getNNodes());
    }

    /**
     * Construct cover from specified clustering.
     *
     * <p>
     * Each node i is assigned to a single cover, as indicated by its cluster.
     * If the cluster of a node is negative, the node is not assigned to any
     * cover.
     * </p>
     *
     * @param clustering Clustering
     * @param nNodes Number of nodes
     */
    public Cover(Clustering clustering, int nNodes)
    {
        this.coverMap = new Int2DoubleOpenHashMap[nNodes];
        for (int i = 0; i < nNodes; i++)
        {
            this.coverMap[i] = new Int2DoubleOpenHashMap();
        }

        for (int i = 0; i < clustering.getNNodes(); i++)
        {
            int cluster = clustering.getCluster(i);
            if (cluster >= 0) // If cluster < 0 it will be "missing".
            {
                this.coverMap[i].put(clustering.getCluster(i), 1);
                if (cluster + 1 > this.nCovers)
                {
                    this.nCovers = cluster + 1;
                }
            }
            coverMap[i].trim();
        }
        this.nNodes = nNodes;
    }

    public Cover clone()
    {
        Cover clonedCovering;
        try
        {
            clonedCovering = (Cover) super.clone();
            clonedCovering.coverMap = coverMap.clone();
            return clonedCovering;
        } catch (CloneNotSupportedException e)
        {
            return null;
        }
    }

    /**
     * Saves Cover object to file.
     *
     * <p>
     * Can be loaded again with load function.
     * </p>
     *
     * @param fileName File to save Cover to
     * @throws IOException Could not write to file
     * @see #load(String filename)
     */
    public void save(String fileName) throws IOException
    {
        ObjectOutputStream objectOutputStream;
        objectOutputStream = new ObjectOutputStream(new FileOutputStream(fileName));
        objectOutputStream.writeObject(this);
        objectOutputStream.close();
    }

    /**
     * Gets number of nodes.
     *
     * @return Number of nodes
     */
    public int getNNodes()
    {
        return nNodes;
    }

    /**
     * Gets number of covers.
     *
     * @return Number of covers
     */
    public int getNCovers()
    {
        return nCovers;
    }

    /**
     * Recalculates the number of covers.
     *
     * <p>
     * This sets the number of covers equal to the minimum cover + 1.
     * </p>
     */
    public void recalcNCovers()
    {
        nCovers = 0;
        for (int i = 0; i < this.nNodes; i++)
        {
            for (int key : this.coverMap[i].keySet())
            {
                if (key + 1 > nCovers)
                {
                    nCovers = key + 1;
                }
            }
        }
    }

    /**
     * Gets the number of covers per node.
     *
     * @return Number of covers per node.
     *
     * @see #getNCovers(int node)
     */
    public int[] getNCoversPerNode()
    {
        int[] nCovers = new int[this.nNodes];
        for (int i = 0; i < this.nNodes; i++)
        {
            nCovers[i] = this.coverMap[i].size();
        }
        return nCovers;
    }

    /**
     * Gets the number of covers of a node.
     *
     * @param node Node
     * @return Number of covers of node
     *
     * @see #getNCoversPerNode()
     */
    public int getNCovers(int node)
    {
        return coverMap[node].size();
    }

    /**
     * Gets the cover weights of node.
     *
     * <p>
     * If node is not inside the specified cover, the cover weight will
     * be 0.0.
     * </p>
     *
     * @param node Node
     * @param cover Cover
     * @return Cover weight
     */
    public double getCoverWeight(int node, int cover)
    {
        return coverMap[node].getOrDefault(cover, 0.0);
    }

    /**
     * Gets sum of node weights for covers for network.
     *
     * <p>
     * The weight of a cover is defined as the sum of all the node weights for
     * all nodes that are assigned to a certain cover.
     * </p>
     *
     * @param network Network
     * @return Cover weights
     *
     * @see #getCoverWeights(double[] nodeWeight)
     */
    public double[] getCoverWeights(Network network)
    {
        return getCoverWeights(network.nodeWeights);
    }

    /**
     * Gets sum of node weights for covers.
     *
     * <p>
     * The weight of a cover is defined as the sum of all the node weights for
     * all nodes that are assigned to a certain cover.
     * </p>
     *
     * @param nodeWeight Node weights
     * @return Cover weights
     *
     * @see #getCoverWeights(Network network)
     */
    public double[] getCoverWeights(double[] nodeWeight)
    {
        double[] coverWeight = new double[nCovers];
        for (int i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[i].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double weight = entry.getDoubleValue();
                coverWeight[cover] += nodeWeight[i] * weight;
            }
        }
        return coverWeight;
    }

    /**
     * Get iterator over all covers for node.
     *
     * <p>
     * This can be used to iterate quickly over all covers using a
     * {@code foreach} structure.
     * </p>
     *
     * @param node Node
     * @return Iterator over all covers for node
     */
    public Int2DoubleMap.FastEntrySet getCoverIterator(int node)
    {
        return this.coverMap[node].int2DoubleEntrySet();
    }

    /**
     * Gets number of nodes for all covers.
     * 
     * @return Number of nodes per cover
     */
    public int[] getNNodesPerCover()
    {
        int i;
        int[] nNodesPerCover;
        nNodesPerCover = new int[nCovers];
        for (i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[i].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                nNodesPerCover[cover]++;
            }
        }
        return nNodesPerCover;
    }

    /**
     * Gets nodes per cover.
     *
     * <p>
     * It so happens that the nodes per cover can also be represented by a
     * cover. Each cover is then considered a node, and each node considered a
     * cover. Essentially, the cover is a sparse matrix representation of the
     * cover weight membership matrix S_{ic}, where S_{ic} =
     * {@code getCoverWeight(i, c)}. The nodes per cover is simply the transpose
     * of S_{ic}.
     * </p>
     *
     * @return Nodes per cover
     */
    public Cover getNodesPerCover()
    {
        int i;
        Int2DoubleOpenHashMap[] nodePerCoverMap = new Int2DoubleOpenHashMap[nCovers];
        for (i = 0; i < nCovers; i++)
        {
            nodePerCoverMap[i] = new Int2DoubleOpenHashMap();
        }

        for (i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[i].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                nodePerCoverMap[cover].put(i, coverWeight);
            }
        }
        Cover cover = new Cover(nodePerCoverMap);
        cover.trim();
        return cover;
    }

    /**
     * Creates subcovers for all clusters in clustering.
     *
     * <p>
     * Generates a cover for each cluster, where each cover
     * only contains the nodes in the specified cluster.
     * </p>
     * 
     * @param clustering Clustering
     * @return Subcovers
     *
     * @see #createSubcover(int[] nodes)
     */
    public Cover[] createSubcovers(Clustering clustering)
    {
        int[][] nodePerCluster;
        Cover[] subcovers;
        nodePerCluster = clustering.getNodesPerCluster();
        subcovers = new Cover[clustering.getNClusters()];
        for (int c = 0; c < clustering.getNClusters(); c++)
        {
            subcovers[c] = this.createSubcover(nodePerCluster[c]);
        }
        return subcovers;
    }

    /**
     * Creates subcover for specified nodes.
     *
     * <p>
     * Generates a cover that only contains the specified nodes.
     * </p>
     *
     * <p>
     * The internal cover map is not cloned, so that the subcover will share the
     * cover map with the overall cover. This is intended behavior, and allows
     * to make changes to subcovers that are automatically reflect in the
     * overall cover.
     * </p>
     *
     * <p>
     * For performance reasons the number of covers is simply set to the same
     * number of covers for the overall cover.
     * </p>
     *
     * @param nodes Nodes
     * @return Subcover
     * 
     * @see #createSubcovers(Clustering clustering)
     * @see #createSubcoverRelabel(int[] nodes, int coverToRelabelToZero)
     */
    public Cover createSubcover(int[] nodes)
    {
        Cover subcover = new Cover(nodes.length, false);
        subcover.coverMap = new Int2DoubleOpenHashMap[nodes.length];
        for (int i = 0; i < nodes.length; i++)
        {
            subcover.coverMap[i] = this.coverMap[nodes[i]]; // Not set clone for memory reasons
        }
        subcover.nCovers = this.nCovers; // Just put it to this.nCovers, that will be the maximum.
        return subcover;
    }

    /**
     * Creates subcover for specified nodes, while relabeling them.
     *
     * <p>
     * The relabeling is done to ensure all covers are in the range of
     * {@code 0, ..., nCovers} where {@code nCovers} is the minimum number of
     * covers required to represent all covers. Moreover, the cover that will be
     * relabeled as zero is indicated by the {@code coverToRelableToZero}. If
     * this cover is not present, nCovers will be one larger than the minimum
     * number of required covers.
     * </p>
     *
     * @param nodes Nodes
     * @param coverToRelabelToZero Cover to relabel to zero
     * @return Subcover
     *
     * @see #createSubcover(int[] nodes)
     */
    public Cover createSubcoverRelabel(int[] nodes, int coverToRelabelToZero)
    {
        Cover subcover = new Cover(nodes.length, true);
        Int2IntOpenHashMap newCoverIds = new Int2IntOpenHashMap();
        newCoverIds.put(coverToRelabelToZero, 0);
        for (int i = 0; i < nodes.length; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[nodes[i]].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                if (!newCoverIds.containsKey(cover))
                {
                    newCoverIds.put(cover, newCoverIds.size());
                }
                subcover.setCoverWeight(i, newCoverIds.get(cover), coverWeight);
            }
        }
        subcover.trim();
        subcover.nCovers = newCoverIds.size();
        return subcover;
    }

    /**
     * Calculates the total overlap of the cover.
     *
     * <p>
     * The total overlap is the sum of the overlap of all nodes.
     * </p>
     * 
     * @return Total overlap
     *
     * @see #calcOverlapSquare(int node)
     * @see #calcOverlapSquare(int node, Cover otherCover, int otherNode)
     */
    public double calcOverlapSquare()
    {
        double overlapSquare = 0.0;
        for (int i = 0; i < this.nNodes; i++)
        {
            overlapSquare += calcOverlapSquare(i);
        }
        return overlapSquare;
    }

    /**
     * Calculates the overlap of the node.
     *
     * <p>
     * The overlap is defined as sum_c S_{ic}S_{ic} over all covers c.
     * </p>
     *
     * @param node Node
     * @return Overlap
     *
     * @see #calcOverlapSquare()
     * @see #calcOverlapSquare(int node, Cover otherCover, int otherNode)
     */
    public double calcOverlapSquare(int node)
    {
        double overlapSquare = 0.0;

        for (Int2DoubleMap.Entry entry : coverMap[node].int2DoubleEntrySet())
        {
            double coverWeight = entry.getDoubleValue();
            overlapSquare += coverWeight * coverWeight;
        }
        return overlapSquare;
    }

    /**
     * Calculates the overlap with another cover.
     *
     * <p>
     * The overlap is defined as sum_c S_{ic} S'_{jc} where S_{ic} is the cover
     * weight of node i for cover c of this cover, and S'_{jc} is the cover
     * weight of the other node j for cover c in the other cover.
     * </p>
     *
     * @param node Node
     * @param otherCover Other cover
     * @param otherNode Node in other cover
     * @return Overlap
     */
    public double calcOverlapSquare(int node, Cover otherCover, int otherNode)
    {
        double overlapSquare = 0.0;

        if (this.getNCovers(node) > otherCover.getNCovers(otherNode))
        {
            for (Int2DoubleMap.Entry entry : this.coverMap[node].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                overlapSquare += coverWeight * otherCover.getCoverWeight(otherNode, cover);
            }
        }
        else
        {
            for (Int2DoubleMap.Entry entry : otherCover.getCoverIterator(otherNode))
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                overlapSquare += coverWeight * this.coverMap[node].get(cover);
            }
        }
        return overlapSquare;
    }

    /**
     * Sets the cover weight for specified node and cover.
     *
     * @param node Node
     * @param cover Cover
     * @param coverWeight Cover weight
     */
    public void setCoverWeight(int node, int cover, double coverWeight)
    {
        this.coverMap[node].put(cover, coverWeight);
    }

    /**
     * Adds weight to cover weight for specified node and cover.
     * 
     * @param node Node
     * @param cover Cover
     * @param coverWeight Cover weight
     */
    public void addToCoverWeight(int node, int cover, double coverWeight)
    {
        this.coverMap[node].addTo(cover, coverWeight);
    }

    /**
     * Resets all covers for node.
     *
     * <p>
     * The node is not assigned to any cover.
     * </p>
     *
     * @param node node
     */
    public void resetCovers(int node)
    {
        this.coverMap[node].clear();
    }

    /**
     * Initializes singleton cover.
     *
     * <p>
     * A singleton cover is such that each node is only in its own cover.
     * </p>
     */
    public void initSingletonCovers()
    {
        initSingletonCoversHelper();
    }

    /**
     * Removes empty covers.
     *
     * <p>
     * Relabels all covers to ensure they are in the range of
     * {@code 0, ..., nCovers} where {@code nCovers} is the minimum number of
     * covers required to represent all non-empty covers.
     * </p>
     *
     */
    public void removeEmptyCovers()
    {
        boolean[] nonEmptyCover;
        int i, j;
        int[] newCover;
        nonEmptyCover = new boolean[nCovers];
        newCover = new int[nCovers];
        for (i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[i].int2DoubleEntrySet())
            {
                nonEmptyCover[entry.getIntKey()] = true;
            }
        }

        i = 0;
        for (j = 0; j < nCovers; j++)
        {
            if (nonEmptyCover[j])
            {
                newCover[j] = i;
                i++;
            }
        }
        nCovers = i;
        Int2DoubleOpenHashMap[] newCoverMap = new Int2DoubleOpenHashMap[this.nNodes];
        for (i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : coverMap[i].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                newCoverMap[i].put(newCover[cover], coverWeight);
            }
            newCoverMap[i].trim();
        }
        this.coverMap = newCoverMap;
    }

    /**
     * Trim internal hashmaps.
     *
     * <p>
     * Trimming the hashmaps makes it as small as possible, thereby reducing the
     * memory footprint.
     * </p>
     *
     * @return Whether trimming was successful
     */
    public boolean trim()
    {
        boolean update = false;
        for (int i = 0; i < this.coverMap.length; i++)
        {
            update |= this.coverMap[i].trim();
        }
        return update;
    }

    /**
     * Create reduced cover based on specified clustering.
     *
     * <p>
     * The reduced cover is created such that the cover weight for cluster c for
     * cover d is the sum of all cover weight for cover d of all node i in
     * cluster c. In other words, S'_{cd} = sum_{i} S_{id} d(s_i, c), where
     * S'_{cd} is the weight of cover d for cluster c in the reduced cover,
     * S_{id} is the weight of cover d for node i and d(s_i, c)=1 if node i
     * is in cluster c and zero otherwise.
     * </p>
     * 
     * @param clustering Clustering
     * @return Reduced cover
     */
    public Cover createReducedCover(Clustering clustering)
    {
        Cover reducedCover = new Cover(clustering.getNClusters());
        for (int i = 0; i < nNodes; i++)
        {
            for (Int2DoubleMap.Entry entry : this.coverMap[i].int2DoubleEntrySet())
            {
                int cover = entry.getIntKey();
                double coverWeight = entry.getDoubleValue();
                reducedCover.addToCoverWeight(clustering.getCluster(i), cover, coverWeight);
                if (cover + 1 > reducedCover.nCovers)
                {
                    reducedCover.nCovers = cover + 1;
                }
            }
        }
        reducedCover.trim();
        return reducedCover;
    }

    private void initSingletonCoversHelper()
    {
        int i;
        this.coverMap = new Int2DoubleOpenHashMap[this.nNodes];
        for (i = 0; i < this.nNodes; i++)
        {
            this.coverMap[i] = new Int2DoubleOpenHashMap();
            this.coverMap[i].put(i, 1);
            this.coverMap[i].trim();
        }
        nCovers = nNodes;
    }
}

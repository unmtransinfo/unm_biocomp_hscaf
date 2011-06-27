package edu.unm.health.biocomp.hscaf;

import chemaxon.struc.MoleculeGraph;
//import chemaxon.util.BondTable;	//API in JChem 5.4
//import chemaxon.core.util.BondTable;	//API change in JChem 5.5

/**
 * The <code>ShortestPath</code> class calculates and stores the length 
 * of shortest paths between any two nodes of a molecular graph&#46; The 
 * implementation is based on the Floyd-Warshall algorithm&#46;
 * One instance of <code>ShortestPath</code> can be reused to compute distances 
 * for several chemical graphs, this serves the more efficient storage space 
 * usage&#46; Internal structures are grown dynamically when a molecular structure 
 * larger than any previous one is given as input for the shortest path 
 * calculation&#46; (Author:  Miklos Vargyas, ChemAxon)
 * <p>
 * <i> Example of typical usage: </i>
 * <pre>
 * ShortestPath sp();   
 * &#46;&#46;&#46;
 * MoleculeGraph mg = get the first molecular graph
 * while ( &#46;&#46;&#46;  there are more molecular graphs &#46;&#46;&#46; ) {
 *   sp&#46;calculate( mg );
 *   for ( int i = 1; i &lt; mg&#46;getAtomCount(); ++i ) {
 *     for ( int j = 0; j &lt; i; ++j ) {
 *       System&#46;out&#46;println( "The length of the shortest path between nodes " +
 *                           i + " and " + j " is " sp&#46;minDist( i, j ) );
 *     }
 *   }
 *   mg = get next molecular graph
 * }
 * </pre>
 *
 * In its current implementation the edges of the graph are unweighted (that is 
 * all weights are equal to 1)&#46; This can easily be changed if needed later&#46; 
 *
	<blockquote>
 	<b>NOTE:</b> This code, downloaded from the ChemAxon forum, has been slightly revised
	to comply with the newer JChem 5&#46;4&#46;0 API&#46;  CGraph has been replaced by MoleculeGraph&#46;
	However, although MoleculeGraph&#46;getBondTable() is recommended, this
	chemaxon.util.BondTable is not documented and it is not clear how to
	access the values, so MoleculeGraph&#46;getBtab() is still used&#46;
	<br>
	<i>(Jeremy Yang, 25 Feb 2011)</i>
	</blockquote>
 *
 *
 * @author  Miklos Vargyas
 * @since JChem   2.0.0
 */
public class ShortestPath {

    /** number of nodes in the current graph */
    protected int nNodes = 0;   
    
    /** maximum number of nodes used (size of storage allocated) */
    protected int maxNodes = 0; 
    
    // temporaries used by the Floyd-Warshall
    
    /** 
     * for each node stores the index of the index predecessor 
     * node on the shortest path found
     */
    private int[][] pred = null;
    
    /** for each node it stores the length of the shortest path found */
    private int[][] d = null;
    
    /** as returned by MoleculeGraph.getBtab() */
    protected int[][] bTab = null;
    /** as returned by MoleculeGraph.getBondTable() */
    //protected BondTable bTab = null; // Use this when we learn how...
                             
    /** 
     * Creates a new ShortestPath object. 
     */
    public ShortestPath() 
    {
        maxNodes = 0;
        nNodes = 0;
    }
    
    /** 
     * Calculates all minimum path lengths and stores them.
     * @param mg The chemical graph in which shortest paths have to be determined.
     */
    public void calculate( MoleculeGraph mg ) 
    {
        nNodes = mg.getAtomCount();
        bTab = mg.getBtab(); // Deprecated...
        //bTab = mg.getBondTable(); // Use this when we learn how...
        alloc();
        init();
        
        for ( int k = 0 ; k < nNodes; k++ ) {
            for ( int j = 0; j < nNodes; j++ ) {
                for ( int i = 0; i < nNodes; i++ ) {
                    int s = d[ i ][ k ] + d[ k ][ j ];
                    if ( s < d[ i ][ j ] ) {
                        d[ i ][ j ] = s;
                        pred[ i ][ j ] = pred[ k ][ j ];
                    }
                }
            }
        }
    }
    
    /**
     * Gets the length of the shortest path between two nodes of the graph.
     * @param  a  index of the first node
     * @param  b  index of the second node
     * @return    The minimum distance (the length of the shortest path) between
     *            the two given nodes.
     */
    public int minDist( int a, int b ) 
    {
        // assert 0 <= a && a < nNodes && 0 <= b && b < nNodes;
	
        return d[ a ][ b ];
    }
    
    /**
     * Gets all nodes on the path between the two givennodes.
     * @param   a   starting node
     * @param   b   end node
     * @param   nodes   output parameter, node indices on the path from a to b
     *                  It is the callers responsibility to provide a 
     *                  sufficiently large array.
     * @return  path length, number of node indices stored in nodes[]
     */
    public int getPath( int a, int b, int[] nodes ) {
        int len = 0;
        while ( pred[ a ][ b ] != -1 && a != b ) {
            nodes[ len++ ] = b;
            b = pred[ a ][ b ];
        }
        nodes[ len++ ] = a;
        return len;
    }
    
    /**
     * Initialises internal structures.
     */ 
    private void init()
    {
        for ( int i = 0; i < nNodes; ++i ) {
            for ( int j = i + 1; j < nNodes; ++j ) {
                if ( bTab[ i ][ j ] != -1 ) {
                    d[ i ][ j ] = 1;
                    d[ j ][ i ] = 1;
                    pred[ i ][ j ] = i;
                    pred[ j ][ i ] = j;
                }
                else {
                    d[ i ][ j ] = nNodes;
                    d[ j ][ i ] = nNodes;
                    pred[ i ][ j ] = -1;
                    pred[ j ][ i ] = -1;
                }
            }
        }
    }
    

    /**
     * Allocates/reallocates internal storage for the current <code>maxNodes</code>
     * value.
     */
    private void alloc()
    {
        if ( nNodes > maxNodes ) {
            maxNodes = nNodes;
            pred = new int[ nNodes ][ nNodes ];
            d = new int[ nNodes ][ nNodes ];
        }
    }
    
}

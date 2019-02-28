/*-----------------------------------------------------------------------

  File        : VertexCovering.h

  Description : Vertex covering by paths on trees

  This is based on the paper "Two fixed-parameter algorithms for
  Vertex Covering by Paths on Trees" by
  Jiong Guo, Rolf Niedermeier and Johannes Uhlmann, in Information
  Processing Letters 106(2008) 81-86

  Copyright  : Francois-Xavier Dupé - http://www.greyc.ensicaen.fr/~fdupe/
               	       
  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.

  ---------------------------------------------------------------------*/

#if !defined(VERTEXCOVERING_H)
#define VERTEXCOVERING_H

#include <pandore.h>
#include <map>
#include <deque>
#include "BagOfTrails.h"
#include "CImg.h"

class VertexCovering
{
private:

  /**
   * Structure needed for the two dynamic algorithms for vertex covering
   */
  struct NodeInfo
  {
    /**
     * vpath(v) contains the subset collection for the node v
     */
    std::deque<trail> vpath;

    /**
     * wpath(v) contains the weight of the paths contained in vpath(v)
     */
    std::deque<double> wpath;

    /**
     * epath(e) for an edge a contains all paths covering e
     */
    std::deque<trail> epath;

    /**
     * The dynamic vector
     */
    double * dynaTab;
  };

  /**
   * Structure needed for the two dynamic algorithms for edge covering
   */
  struct EdgeInfo
  {
    /**
     * The parent node of the edge
     */
    int parent;
    
    /**
     * The child node of the edge
     */
    int child;

    /**
     * True if the edge is the root edge
     */
    bool root;

    /**
     * wpath(e) contains the weight of the paths contained in epath(e)
     */
    std::deque<double> wpath;

    /**
     * epath(e) for an edge a contains all paths covering e
     */
    std::deque<trail> epath;

    /**
     * The dynamic vector
     */
    double * dynaTab;
  };

  /**
   * The current graph
   */
  pandore::Collection * col;
  
  /**
   * The current bag of paths
   */
  std::deque<trail> bag;

  /**
   * The solution (kept it in memory)
   */
  std::deque<trail> result;
  
  /**
   * Get a suffix trip of the graph (using the first non-null node as root)
   * @return a suffix representation of the tree
   */
  std::deque<int> getVertexSuffixTree ( void );

  /**
   * Recursive version for suffix computation (vertex version)
   * @param grp    the graph (here a tree)
   * @param tab    the current trip
   * @param node   the current node
   * @param visit  the visit table
   */
  void getVertexSuffixTree ( pandore::Graph3d * grp, std::deque<int> & tab, int node, cimg_library::CImg<bool> & visit );

  /**
   * Compute the vpath sets for each node
   * @param grp    the graph
   * @param bag    the input bag of path
   * @param suffix the suffix representation of the tree
   * @param sNode  the array of nodes' information
   */
  void computeVpath ( pandore::Graph3d * grp, const std::deque<trail> & bag,
		      const std::deque<int> & suffix, std::map<int,NodeInfo> & sNode );

  /**
   * Compute the union of two sets of paths
   * @param s1 the first set
   * @param s2 the second set
   * @return the strict union of the two sets
   */
  std::deque<trail> unionset ( const std::deque<trail> & s1, const std::deque<trail> & s2 );

  /**
   * Compute the intersection of two sets of paths
   * @param s1 the first set
   * @param s2 the second set
   * @return the intersection of the two sets
   */
  std::deque<trail> intersection ( const std::deque<trail> & s1, const std::deque<trail> & s2 );

  /**
   * Compute the assymetric difference between two sets of paths
   * @param s1 the first set
   * @param s2 the second set
   * @return the difference between the two sets
   */
  std::deque<trail> difference ( const std::deque<trail> & s1, const std::deque<trail> & s2 );

  /**
   * Compute the maximal weight vertex covering by paths of the tree
   * @param suffix  the suffix representation of the tree
   * @param sNode   the information about each node
   * @return the set of paths which covers the tree
   */
  std::deque<trail> pathVertexCovering ( const std::deque<int> & suffix, std::map<int,NodeInfo> & sNode );

  /**
   * Compute the vertex solution using the dynamic table results
   * @param node   the current node
   * @param sNode  the dynamic table
   * @param res    the optimal set of paths
   * @param visit  the visit table
   */
  void computeVertexSolution ( int node, std::map<int,NodeInfo> & sNode, std::deque<trail> & res,
			       cimg_library::CImg<bool> & visit );

public:

  /**
   * The constructor
   * @param graph  the graph to cover
   * @param path   the bag of path
   */
  VertexCovering ( pandore::Collection * col, const std::deque<trail> & bag );

  /**
   * Compute covering
   * @param edge   do an edge covering instead of a vertex covering
   * @param lambde the parameter setting the equilibrium between cardinalty and weight
   * @return the covering set
   */
  std::deque<trail> computeCovering ( bool edge, double lambda );
};

#endif

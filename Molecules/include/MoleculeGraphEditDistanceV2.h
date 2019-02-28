/**
 * @file MoleculeGraphEditDistanceV2.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __MOLECULE_GRAPH_EDIT_DISTANCE_V2_H__
#define __MOLECULE_GRAPH_EDIT_DISTANCE_V2_H__

#include <deque>
#include <string>
#include <map>

#include "pandore.h"
#include "CImg.h"
#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "GraphEditDistance.h"
#include "Dataset.h"

/**
 * @brief The graph edit distance between two molecules.
 *
 * This is the specialization of the graph edit distance for molecules.
 *
 * Since a node substitution/addition/deletion can have important consequences
 * on the activity of a molecule, these operations must be strongly penalized.
 *
 * The operations costs are :
 *		- Node substitution : 0 if same label (atom symbol), 1 otherwise.
 *				+ 2 for each extra neighbour
 *		- Node addition/deletion : 2 + number of neighbors of the node.
 */


class MoleculeGraphEditDistanceV2 : public GraphEditDistance
{
public:
  

  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  
  /**
   * The node substitution cost.
   *
   * @param node1 The first node.
   * @param node2 The second node.
   *
   * @return The cost of the substitution node1 -> node2.
   */
	
  int nodeSubstitutionCost(pandore::Collection* node1, pandore::Collection* node2)
  {
    return !(node1->GETVALUE("atom", pandore::Char) == node2->GETVALUE("atom", pandore::Char));
  }
	
  /**
   * The adjacency substitution cost.
   *
   * @param node1 The first node.
   * @param node2 The second node.
   *
   * @return The cost of the substitution of node1 -> node2.
   */
	
  int adjacencySubstitutionCost(pandore::Collection* node1, pandore::Collection* node2)
  {
    return 0;
  }
	
  /**
   * The node addition/deletion cost.
   *
   * @param node1 The node to insert/delete.
   * @param e The node's neighbours.
   *
   * @return The cost of the node's insertion/deletion.
   */
	
  int extraNodeCost (pandore::Collection* node, pandore::GEdge* e);
	
  /**
   * The adjacency addition/deletion cost.
   *
   * @param node1 The node to insert/delete.
   *
   * @return The cost of the node's insertion/deletion.
   */
	
  int extraNeighbourCost (pandore::Collection* node)
  {
    return 2;
  }

       
};

#endif // __MOLECULE_GRAPH_EDIT_DISTANCE_V2_H__

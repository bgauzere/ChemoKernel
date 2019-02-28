/**
 * @file GraphEditDistance.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __GRAPH_EDIT_DISTANCE_H__
#define __GRAPH_EDIT_DISTANCE_H__

#include "EditDistance.h"
#include <pandore.h>

/**
 * @brief Approximation of the graph edit distance
 *
 * This class implements the approximation of the edit
 * distance between two graphs, using the Munkres algorithm.
 */

class GraphEditDistance : public EditDistance
{
public:
	
	/**
	 * Returns the approximated edit distance between the two graphs c1 and c2.
	 *
	 * @param c1 The first graph.
	 * @param c2 The second graph.
	 *
	 * @return The approximated edit distance between c1 and c2.
	 */
  
         double operator() (pandore::Collection* c1, pandore::Collection* c2);
	
	/**
	 * The node substitution cost.
	 *
	 * @param node1 The first node.
	 * @param node2 The second node.
	 *
	 * @return The cost of the substitution node1 -> node2.
	 */
	
	virtual int nodeSubstitutionCost(pandore::Collection* node1, pandore::Collection* node2) = 0;
	
	/**
	 * The adjacency substitution cost.
	 *
	 * @param node1 The first node.
	 * @param node2 The second node.
	 *
	 * @return The cost of the substitution of node1 -> node2.
	 */	
	
	virtual int adjacencySubstitutionCost(pandore::Collection* node1, pandore::Collection* node2) = 0;
	
	/**
	 * The node addition/deletion cost.
	 *
	 * @param node1 The node to insert/delete.
	 * @param e The node's neighbours.
	 *
	 * @return The cost of the node's insertion/deletion.
	 */
	
	virtual int extraNodeCost (pandore::Collection* node, pandore::GEdge* e) = 0;
	
	/**
	 * The adjacency addition/deletion cost.
	 *
	 * @param node1 The node to insert/delete.
	 *
	 * @return The cost of the node's insertion/deletion.
	 */
	
	virtual int extraNeighbourCost (pandore::Collection* node) = 0;
};

#endif // __GRAPH_EDIT_DISTANCE_H__

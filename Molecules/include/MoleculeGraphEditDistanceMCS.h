/**
 * @file MoleculeGraphEditDistanceMCS.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __MOLECULE_GRAPH_EDIT_DISTANCE_MCS_H__
#define __MOLECULE_GRAPH_EDIT_DISTANCE_MCS_H__

#include <deque>
#include <string>
#include <map>

#include "pandore.h"
#include "CImg.h"
#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "GraphEditDistance.h"
#include "Dataset.h"

class MoleculeGraphEditDistanceMCS : public GraphEditDistance
{
  double c_s;//cout de substitution
  double c_i;//cout de modif structurelle

public:
  
  MoleculeGraphEditDistanceMCS(double c_s,double c_i);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  int nodeSubstitutionCost(pandore::Collection* node1, pandore::Collection* node2)	{
    return 0;
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
	
  int extraNodeCost (pandore::Collection* node, pandore::GEdge* e)	{
    return 0;
  }

	
  /**
   * The adjacency addition/deletion cost.
   *
   * @param node1 The node to insert/delete.
   *
   * @return The cost of the node's insertion/deletion.
   */
	
	
  int extraNeighbourCost (pandore::Collection* node)
  {
    return 0;
  }

  

       
};

#endif // __MOLECULE_GRAPH_EDIT_DISTANCE_MCS_H__

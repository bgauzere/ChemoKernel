/*
 * @file MoleculeGraphEditDistanceV2.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "MoleculeGraphEditDistanceV2.h"

using namespace pandore;
using namespace cimg_library;


double MoleculeGraphEditDistanceV2::operator()(Collection* c1, Collection* c2)
{
  double lambda = 1.0;
  Collection** nodes1 = c1->GETPARRAY("nodes", Collection);
  Collection** nodes2 = c2->GETPARRAY("nodes", Collection);
  
  Graph3d* graph1 = c1->GETPOBJECT("graph", Graph3d);
  Graph3d* graph2 = c2->GETPOBJECT("graph", Graph3d);
  
  // Computing the node cost matrix
  unsigned int n1 = c1->GETPARRAYSIZE("nodes", Collection);
  unsigned int n2 = c2->GETPARRAYSIZE("nodes", Collection);
  unsigned int n = n1;	
  if (n2 > n)
    n = n2;

  CImg<double> node_cost (n, n, 1, 1, 0);
  for (unsigned int i=0; i<n1; ++i)
    {
      for (unsigned int j=0; j<n2; ++j)
  	{
  	  bool central_node_similar = !(nodes1[i]->GETVALUE("atom", pandore::Char) == 
  					nodes2[j]->GETVALUE("atom", pandore::Char));
  	  //Number of neighbours
  	  unsigned int p1 = 0; // The number of neighbours of the node i
  	  for (GEdge* e1 = (*graph1)[i]->Neighbours(); e1!=NULL; e1=e1->Next())
  	    ++p1;
  	  unsigned int p2 = 0; // The number of neighbours of the node j
  	  for (GEdge* e2 = (*graph2)[j]->Neighbours(); e2!=NULL; e2=e2->Next())
  	    ++p2;
  	  unsigned int p = (p1>p2)?p1:p2;
						
  	  // The edge cost matrix is computed
  	  cimg_library::CImg<double> edge_cost (p, p, 1, 1, 0);
  	  GEdge* e1 = (*graph1)[i]->Neighbours();
  	  for (unsigned int k=0; k<p; ++k)
  	    {
  	      GEdge* e2 = (*graph2)[j]->Neighbours();
  	      for (unsigned int l=0; l<p; ++l)
  		{
  		  bool neighbours_are_similar = false;
  		  if ((e1 != NULL) && (e2 != NULL))
  		    neighbours_are_similar = !(nodes1[e1->Node()]->GETVALUE("atom", pandore::Char) == 
  					       nodes2[e2->Node()]->GETVALUE("atom", pandore::Char));
  		  edge_cost(k,l) +=  ((lambda*central_node_similar) + neighbours_are_similar);
		  if (e2 != NULL)
  		    e2 = e2->Next();
  		}
  	      if (e1 != NULL)
  		e1 = e1->Next();
  	    }
			
  	  // We apply Munkres algorithm on the adjacency cost matrix
  	  Alignement* align = applyMunkres(edge_cost);
  	  node_cost(i,j) = align->cost;
  	  delete align;
  	}
      //Appariemment 
      for (unsigned int j=n2; j<n; ++j)
  	node_cost(i,j) += extraNodeCost(nodes1[i], (*graph1)[i]->Neighbours());
    }
	
  for (unsigned int i=n1; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      node_cost(i,j) += extraNodeCost(nodes2[j], (*graph2)[j]->Neighbours());
  
  // We apply Munkres algorithm on the cost matrix
  Alignement* align = applyMunkres(node_cost);
  double cost = align->cost;
  delete align;
	
  return cost;
}

int MoleculeGraphEditDistanceV2::extraNodeCost (pandore::Collection* node, pandore::GEdge* e)
{
  int p=0;
  while (e != NULL)
    {
      ++p;
      e = e->Next();
    }
  return 2*p;
}

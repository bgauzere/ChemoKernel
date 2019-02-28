/*
 * @file MoleculeGraphEditDistanceMCS.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include <cfloat>
#include "MoleculeGraphEditDistanceMCS.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

MoleculeGraphEditDistanceMCS::MoleculeGraphEditDistanceMCS(double c_s,double c_i){
  this->c_s=c_s;
  this->c_i=c_i;

}

double MoleculeGraphEditDistanceMCS::operator()(Collection* c1, Collection* c2)
{
  Collection** nodes1 = c1->GETPARRAY("nodes", Collection);
  Collection** nodes2 = c2->GETPARRAY("nodes", Collection);
  char* edges1 = c1->GETARRAY("edges", Char);
  char* edges2 = c2->GETARRAY("edges", Char);

  Graph3d* graph1 = c1->GETPOBJECT("graph", Graph3d);
  Graph3d* graph2 = c2->GETPOBJECT("graph", Graph3d);
  
  // Computing the node cost matrix
  unsigned int n1 = c1->GETPARRAYSIZE("nodes", Collection);
  unsigned int n2 = c2->GETPARRAYSIZE("nodes", Collection);
  unsigned int n = n1;	
  if (n2 > n)
    n = n2;

  int * nb_neighbours1 = new int[n1];
  int * nb_neighbours2 = new int[n2];

  CImg<double> node_cost (n, n, 1, 1, 0);
  for (unsigned int i=0; i<n1; ++i)
    {
      //Number of neighbours
      int p1 = 0; // The number of neighbours of the node i
      for (GEdge* e1 = (*graph1)[i]->Neighbours(); e1!=NULL; e1=e1->Next())
	++p1;
      nb_neighbours1[i] = p1;
      for (unsigned int j=0; j<n2; ++j)
  	{
  	  int p2 = 0; // The number of neighbours of the node j
  	  for (GEdge* e2 = (*graph2)[j]->Neighbours(); e2!=NULL; e2=e2->Next())
  	    ++p2;
	  nb_neighbours2[j] = p2;
  	  int p = (p1<p2)?p1:p2;
	  int diff_nb_neighbours = abs(p1-p2);
  	  // The edge cost matrix is computed
  	  cimg_library::CImg<double> edge_cost (p, p, 1, 1, 0);
  	  GEdge* e1 = (*graph1)[i]->Neighbours();
  	  for (unsigned int k=0; k<p; ++k)
  	    {
  	      GEdge* e2 = (*graph2)[j]->Neighbours();
  	      for (unsigned int l=0; l<p; ++l)
  		{
		  edge_cost(k,l) =  !(nodes1[e1->Node()]->GETVALUE("atom", pandore::Char) == 
				      nodes2[e2->Node()]->GETVALUE("atom", pandore::Char));
		  edge_cost(k,l) += !(edges1[e1->Item()] == edges2[e2->Item()]);
		  
		  if (e2 != NULL)
  		    e2 = e2->Next();
  		}
  	      if (e1 != NULL)
  		e1 = e1->Next();
  	    }
	  // We apply Munkres algorithm on the adjacency cost matrix
  	  Alignement* align = applyMunkres(edge_cost);
	  bool central_node_similar = !(nodes1[i]->GETVALUE("atom", pandore::Char) == 
  					nodes2[j]->GETVALUE("atom", pandore::Char));
  	  node_cost(i,j) = (align->cost + central_node_similar)*c_s + (diff_nb_neighbours*c_i*2);
  	  delete align;
  	}
      //Appariemment 
      for (unsigned int j=n2; j<n; ++j)
	node_cost(i,j) = c_i + c_i*p1 + c_i*p1;
    }
  for (unsigned int i=n1; i<n; ++i) //n1 < n2
    for (unsigned int j=0; j<n; ++j){
      node_cost(i,j) = c_i + (c_i*nb_neighbours2[j]*2);// On associe epsilon de g1 avec chaque noeud de g2
    }
  // We apply Munkres algorithm on the cost matrix
  Alignement* align = applyMunkres(node_cost);
  double cost = align->cost;
  cost = 0.0;
  //Computing edit path cost
  //Nodes
  int * appariemmentG1G2 = new int[n];
  int * appariemmentG2G1 = new int[n];
  for (unsigned int i=0; i<n; ++i){
    int x = align->alX[i];
    int f_x = align->alY[i];
    appariemmentG1G2[x] = f_x;
    appariemmentG2G1[f_x] = x;
    if((x<n1)&&(f_x>= n2))  // Appariemment sur Epsilon
      cost += c_i;
    else if ((x>=n1) && (f_x < n2))  //Noeud de g2 apparié sur epsilon
      cost += c_i;
    else
      cost += !(nodes1[x]->GETVALUE("atom", pandore::Char) == 
  		nodes2[f_x]->GETVALUE("atom", pandore::Char))*c_s;    
  }
  //Edges
  for (unsigned int i=0; i<n1; ++i){
    if(appariemmentG1G2[i] < n2){//Apparié sur un noeud de G2
      for (GEdge* e1 = (*graph1)[i]->Neighbours(); e1!=NULL; e1=e1->Next()){
	if(appariemmentG1G2[e1->Node()] < n2){
	  bool exist = false;
	  for (GEdge* e2 = (*graph2)[appariemmentG1G2[i]]->Neighbours(); e2!=NULL; e2=e2->Next()){
	    //Parcours des edges du noeud apparié
	    if(appariemmentG1G2[e1->Node()] == e2->Node() ){ //Edge correspondant trouvé
	      cost += (!(edges2[e2->Item()] == edges1[e1->Item()]))*c_s;
	      exist = true;
	    }
	  }if(!exist)//pas d edge correspondant
	     cost += c_i;
	}
      }
    }
    else{//Appariemment de i de G2 sur epsilon
      cost += (nb_neighbours1[i] - 1) * c_i;
    }
  }
  
  for (unsigned int i=0; i<n2; ++i){
    if(appariemmentG2G1[i] < n1){//Appariemment existant
      for (GEdge* e1 = (*graph2)[i]->Neighbours(); e1!=NULL; e1=e1->Next()){
	if(appariemmentG2G1[e1->Node()] < n1){
	  bool exist = false;
	  for (GEdge* e2 = (*graph1)[appariemmentG2G1[i]]->Neighbours(); e2!=NULL; e2=e2->Next()){
	    if(appariemmentG2G1[e1->Node()] == e2->Node()){
	      //si existe pas besoin de le compter, deja inclus dans G1->G2
	      exist = true;
	    }
	  }
	  if(!exist)
	    cost += c_i;
	}
      }
    }
    else{//Appariemment de i de G2 sur epsilon
      cost += (nb_neighbours2[i] - 1) * c_i;
    }
  }

  delete [] appariemmentG1G2;
  delete [] appariemmentG2G1;
  delete [] nb_neighbours1;
  delete [] nb_neighbours2;
  
  delete align;
  return cost;

}


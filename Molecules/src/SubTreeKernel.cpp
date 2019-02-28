/*
 * @file SubTreeKernel.cpp
 *
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *
 * @version     0.0.1 - Jeu 26 Sep 2013
 * 
 *  
 * Implement the tree pattern kernel [Mahé 2009]
 * and its extension to chiral molecules [Brown 2010]
 *
 */

#include "SubTreeKernel.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include "MoleculeGraph.h"
#include "utils.h"
#include "string_utils.h"
#include "separator.h"

using namespace std;
using namespace pandore;

double SubTreeKernel::computeSubTreeKernel(Graph3d& g1,vector<string> nodes1,vector<string> edges1,
					   Graph3d& g2,vector<string> nodes2,vector<string> edges2,
					   int depthMax,double lambda,bool untilFlag,bool branchFlag)
{
  
  // 1 - create matrices of k_1(r,s)...k_h(r,s)
  // --> subMatrices[0],...,subMatrices[h-1]
  double *** subMatrices;
  subMatrices = new double**[depthMax];
  for(int i = 0; i < depthMax ; i++) 
    {
      subMatrices[i]= new double * [g1.Size()];
      for(int j = 0 ; j < g1.Size() ; j++)
	{
	  subMatrices[i][j] = new double [g2.Size()];
	  for(int k = 0 ; k < g2.Size() ; k++)
	    subMatrices[i][j][k]=0.0;      
	}
    }
	
  // 2 - initialize subMatrices[0]
  vector<int> toupdate1,toupdate2;
  for(int i = 0 ; i <  g1.Size() ; i++)
    for(int j = 0 ; j < g2.Size() ; j++)
      if(nodes1[i].compare(nodes2[j])==0)
	{
	  if(branchFlag)
	    subMatrices[0][i][j] = 1;
	  else
	    subMatrices[0][i][j] = lambda;
	  toupdate1.push_back(i);
	  toupdate2.push_back(j);  
	}
	
  // 3 - compute subMatrices[1]...subMatrices[depthMax]
  double update, temp, ker;	
  for(int depth = 1 ; depth < depthMax ; depth++)				
    for (int istep=0 ; istep < toupdate1.size() ; istep++) 
      {
	int atom1 = toupdate1[istep];
	int atom2 = toupdate2[istep];
	vector<double> tempKernel;
			
	// Detect common neighbor labels

	vector< pair<int,int> > matchableNeighbours;

	for(GEdge * e1 = g1[atom1]->Neighbours();e1!=NULL;e1=e1->Next())
	  for(GEdge * e2 = g2[atom2]->Neighbours();e2!=NULL;e2=e2->Next())
	    if(nodes1[e1->Node()].compare(nodes2[e2->Node()])==0)
	      if(edges1[e1->Item()].compare(edges2[e2->Item()])==0)
		matchableNeighbours.push_back(make_pair(e1->Node(),e2->Node()));
	
	// int maxTuplesSize = min(g1[atom1]->Neighbours(),g2[atom2]->Neighbours());
	// for(int tupleSize=1;tupleSize<=maxTuplesSize;tupleSize++)
	//   {
	    


	//   }

      }
// 		  // A matching pair has been found
// 		  int natoms1 = (*edges1)[atom1][n1].i.size();
// 		  int natoms2 = (*edges2)[atom2][n2].i.size();
// 		  int maxCard = min( natoms1, natoms2); 
// 		  temp = 0.0;
// 	  for(int d = 1; d <= maxCard; d++){
// 	    for(int l = 0 ; l < tuples[natoms1-1][d-1].size() ; l++){
// 	      for(int m = 0 ; m < tuples[natoms2-1][d-1].size() ; m++){
// 		update = 1.0;
// 		for (int toto = 0 ; toto<d ; toto++) {
// 		  update = update * subMatrices[depth-1][ (*edges1)[atom1][n1].i[tuples[natoms1-1][d-1][l][toto] ] ][ (*edges2)[atom2][n2].i[tuples[natoms2-1][d-1][m][toto] ] ];
// 		}
							
// 		if(branchFlag){
// 		  temp += update * pow(lambda,d-1);
// 		}
// 		else{
// 		  temp += update;
// 		}
// 	      }//loop m
// 	    } //loop l
// 	  }//loop d
// 	  tempKernel.push_back(temp);
				
// 	  n1++;
// 	  n2++;
// 	} // end of while( (n1<... loop
			
// 	// Compute kernel update from tempkernel
// 	if (tempKernel.size() == 0)
// 	  // No common neighbour has been found
// 	  ker = 0.0;
// 	else {
// 	  ker = tempKernel[0] + 1;
// 	  for (int n = 1 ; n < tempKernel.size() ; n++){
// 	    if(branchFlag){
// 	      ker = ker * (lambda*tempKernel[n] + 1) + tempKernel[n]*(1 - lambda);
// 	    }
// 	    else{
// 	      ker = ker * (tempKernel[n] + 1);
// 	    }
// 	  }
// 	  ker = ker - 1;
// 	}
			
// 	// update DP matrix
// 	subMatrices[depth][atom1][atom2] = ker;
// 	if(untilFlag){
// 	  subMatrices[depth][atom1][atom2] = 1.0 + subMatrices[depth][atom1][atom2];
// 	}
// 	if(!branchFlag){
// 	  subMatrices[depth][atom1][atom2] = lambda * subMatrices[depth][atom1][atom2];
// 	}
			
//       } //loop i
		
// } //loop depth
				
	
	
	
// 	// 4 - finally : sum subMatrices[depthMax-1] to get the kernel
// 	double kernel = 0.0;
// 	for(int i = 0 ; i < mol1->numAtoms() ; i++){
// 		for(int j = 0 ; j < mol2->numAtoms() ; j++){
// 			if( !(noTotters && (mol1->getAtomByIndex(i)->getKashimaPS() * mol2->getAtomByIndex(j)->getKashimaPS() == 0.0) ) ){ 
// 				kernel += subMatrices[depthMax-1][i][j];
// 			}
// 		}
// 	}
	
// 	if(!branchFlag){
// 		kernel = kernel / pow(lambda,depthMax);
// 	}
	
	
	
// 	return kernel;
  return 0;
}


void SubTreeKernel::noTottersTransform(Graph3d& gToTransform,vector<string> nodes,vector<string> edges,
				       Graph3d& gRes,vector<string> & nodesRes,vector<string> & edgesRes,
				       vector<bool> & trueNode)
{

  vector< int > coord1;
  vector< int > coord2;
  vector< string > edgeLabelTMP;

  int nbAtomNew = 0;

  for(int i = 0 ; i <  gToTransform.Size() ; i++)
    {
      nbAtomNew++;
      for(GEdge * e =gToTransform[i]->Neighbours();e!=NULL;e=e->Next())
	nbAtomNew++;
    }

  gRes.New(nbAtomNew, 0, 0, 0);

  int indiceAtomNewGraph = 0;
    
  // 1st step : create new atoms : 1 per original atom + 1 per original bond
  // for each atom in the moldeule --> create an atom
  for(int i = 0 ; i <  gToTransform.Size() ; i++)
    {
      gRes.Add(indiceAtomNewGraph,indiceAtomNewGraph);
      trueNode.push_back(true);
      coord1.push_back(-1);//coord of the 1st atom of the bond (=-1 in the case of an atom)
      coord2.push_back(i); // coord of the second atom of the bond
      edgeLabelTMP.push_back("");  // bond label
      nodesRes.push_back(nodes[i]);
      indiceAtomNewGraph++;

      // for all the atom's bonds --> create a node
      for(GEdge * e =gToTransform[i]->Neighbours();e!=NULL;e=e->Next())
	{
	  gRes.Add(indiceAtomNewGraph,indiceAtomNewGraph);
	  trueNode.push_back(false);
	  coord1.push_back(i);//coord of the 1st atom of the bond
	  coord2.push_back(e->Node()); // coord of the second atom of the bond
	  edgeLabelTMP.push_back(edges[e->Item()]);  // bond label
	  nodesRes.push_back(nodes[e->Node()]);
	  indiceAtomNewGraph++;
	}
    }

  // 2nd step : create the bonds from the coordinates vector
   int nbLinks = 0;

   for(unsigned int i = 0 ; i < coord1.size() ; i++)
       for(unsigned int j = 0 ; j < coord1.size() ; j++)
	 if(coord2[i] == coord1[j] && coord2[j] != coord1[i])
	   {
	     gRes.Link(i,j,nbLinks,1.0f,false);
	     edgesRes.push_back(edgeLabelTMP[j]);
	     nbLinks++;
	   }
}

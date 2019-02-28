/*
 * @file MoleculeGraphUnlabeledSpectrum.cpp
 * @author Benoit Gauzere <<benoit.gauzere@ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb  2 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Implementation article GBR 2011
 * Unlabeled Treelets, GBR 2011
 *
 */

#include <assert.h>
#include "MoleculeGraph.h"

using namespace std;
using namespace pandore;

int MoleculeGraph::nbPathsSizeK_recursive(GNode<Point3d> * n_start, int k, 
					  GNode<Point3d> * n_parent)
{
  GEdge * v = n_start->Neighbours();
  int nb_paths = 0;
  if(k==1)
    return 1;
  while(v != NULL){
    if((n_parent == NULL) || (v->Node() != n_parent->Item()))
      nb_paths += nbPathsSizeK_recursive((*_graph)[v->Node()]
					 , k-1, n_start);
    v = v->Next();
  }
  return nb_paths;
}

int MoleculeGraph::nbPathsSizeK(int k)
{
  /*
   * Possibilité d'améliorer la complexité en incrémentant dgSpectrum dans 
   * nbPathsSizeK_recursive au fur et a mesure de l'avancement ...
   */
  unsigned int n = _graph->Size();
  int nb_paths = 0;
  for(unsigned int i = 0; i < n; i ++){
    nb_paths += nbPathsSizeK_recursive((*_graph)[i], k, NULL);
  }
  return nb_paths;
}


void MoleculeGraph::computeGraphletSpectrum()
{
  //Acyclic, unlabeled and valence  max = 4 (ex:Alkane family)
  //Complexité a revoir
  _sizeSpectrum = SIZE_SPECTRUM;
  _dgSpectrum = new double[_sizeSpectrum];
  memset(_dgSpectrum, 0,_sizeSpectrum*sizeof(double));
  _dgSpectrum[0] = nbPathsSizeK(1);
  /*On divise par 2 pour eviter de compter les chemins symétriques*/
  _dgSpectrum[1] = nbPathsSizeK(2)/2.0;
  _dgSpectrum[2] = nbPathsSizeK(3)/2.0;
  _dgSpectrum[3] = nbPathsSizeK(4)/2.0;
  _dgSpectrum[4] = nbPathsSizeK(5)/2.0;
  _dgSpectrum[5] = nbPathsSizeK(6)/2.0;
  
  unsigned int n = _graph->Size();
  
  for(unsigned int i=0;i<n;i++)
    {
      int nb_neighbours = nbNeighbours(i);
      if(nb_neighbours == 3)
	{
	  // G4 : 3 star centered on i
	  _dgSpectrum[6] ++;
	  int nb_neighbours_with_plus2_neighbours[3] = {0,0,0};
	  int k = 0;
	  for(GEdge * v = (*_graph)[i]->Neighbours();v!=NULL;v=v->Next())
	    {
	      //G6 : 3 star + leaf
	      _dgSpectrum[7] += nbNeighbours(v->Node()) - 1;
	      nb_neighbours_with_plus2_neighbours[k] = nbNeighbours(v->Node()) -1;
	      k++;
	      //G12 : 1 neighbour with at least 3 neighbours
	      if(nbNeighbours(v->Node()) >= 3)
		_dgSpectrum[12] += 0.5; //XXX: symetriques
	      //G10 : Leaf of a leaf of 3-star
	      for(GEdge * v_prime = (*_graph)[v->Node()]->Neighbours();v_prime!=NULL;v_prime=v_prime->Next())
		_dgSpectrum[10] += nbNeighbours(v->Node()) -1;
	    }
	  //G9
	  for(k=0;k<3;k++)
	    {
	      if(nb_neighbours_with_plus2_neighbours[k] != 0){
		_dgSpectrum[9] += nb_neighbours_with_plus2_neighbours[(k+1)%3];
		_dgSpectrum[9] += nb_neighbours_with_plus2_neighbours[(k+2)%3];
	      }
	
	    }
	}
      if(nb_neighbours == 4)
	{
	  _dgSpectrum[8] ++; //G7, 4 star
	  _dgSpectrum[6] += 4; //Les 4 permutations de 3 star dans un graphlet "4-star"
	  //Parcours de chaque extremité du 4 star
	  for(GEdge * v = (*_graph)[i]->Neighbours();v!=NULL;v=v->Next())
	    {
	      int nb_neighbours_extremity = nbNeighbours(v->Node());
	      _dgSpectrum[7] += 3*(nb_neighbours_extremity-1);//Les 3 star + leaf de chaque permutations
	      _dgSpectrum[11] += nb_neighbours_extremity;
	    } 
	}
      if(nb_neighbours == 5)
	{
	  /*Impossible sur les alkanes !! (max conectivité Carbone = 4)
	    _dgSpectrum[13] ++; //G13, 5 star
	  */
	}
    }
  double nb = 0.0;
  for(int i = 0; i < _sizeSpectrum; i++)
    nb += _dgSpectrum[i];
  //Normalization
  // for(int i = 0; i < _sizeSpectrum; i++)
  //   _dgSpectrum[i] /= (double)nb;
  
  _collection->SETARRAY("spectrum", Double, _dgSpectrum, _sizeSpectrum);
  _collection->SETVALUE("nb_graphlets", Double, nb);
  
}


double * MoleculeGraph::getSpectrum(){
  assert(_dgSpectrum != NULL);
  return _dgSpectrum;
}

int MoleculeGraph::getSpectrumSize(){
  assert(_dgSpectrum != NULL);
  return _sizeSpectrum;
}

/*End Unlabeled Treelets*/

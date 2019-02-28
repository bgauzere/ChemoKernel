/*
 * @file MoleculeGraphLabeledSpectrum.cpp
 * @author Benoit Gauzere <<benoit.gauzere@ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb  2 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Labeled Treelet Distribution
 * 
 *
 */

#include <iostream>
#include <map>
#include <sstream>
#include "MoleculeGraph.h"
#include "string_utils.h"
#include <cassert>
using namespace std;
using namespace pandore;

void MoleculeGraph::computeLabeledTreeletSpectrum()
{ 

  //Nodes/Edges tables construction
  vector<string> edges;
  for(int i=0;i<nbBonds;i++){
    stringstream s;
    s << _edges[i];
    edges.push_back(s.str());
  }
  
  vector<string> nodes;
  for(int i=0;i<this->nbAtoms();i++){
    char c = _nodes[i]->GETVALUE("atom", Char);
    stringstream s;
    s << c;
    nodes.push_back(s.str());
  }
  treelet_spectrum ** distrib = TreeletEnumerator::computeTreeletSpectrum((*_graph),nodes,edges);
  this->labeled_spectrum = distrib;
  long cast_labeled_spectrum = (long)(this->labeled_spectrum);
  assert(sizeof(long) == sizeof(this->labeled_spectrum));
  assert(sizeof(long) == sizeof(Long));  
  
  //XXX : this is too much pointer size dependant. May not work on 32 bits
  _collection->SETVALUE("labeled_spectrum", Long, (Long)(cast_labeled_spectrum));
}

void MoleculeGraph::normalizeLabeledSpectrum()
{
  // double nb_graphlets = 0.0;
  // for(int i = 0; i < SIZE_SPECTRUM; i++){
  //   for(treelet_spectrum::iterator it = labeled_spectrum[i]->begin(); it != labeled_spectrum[i]->end();it++)
  //     nb_graphlets += it->second;
  // }
  
  // for(int i = 0; i < SIZE_SPECTRUM; i++){
  //   for(treelet_spectrum::iterator it = labeled_spectrum[i]->begin(); it != labeled_spectrum[i]->end();it++)
  //     {
  // 	it->second /= nb_graphlets;
  //     }
  // }
}

treelet_spectrum ** MoleculeGraph::getLabeledTreeletSpectrum(){
  return labeled_spectrum;
}

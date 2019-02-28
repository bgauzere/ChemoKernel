/*
 * @file MoleculeGraphAugmentedCycle.cpp
 * @author  <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Oct 17 2013
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include <map>
#include <sstream>
#include <algorithm>
#include <cstring>
#include "MoleculeGraph.h"
#include "string_utils.h"
#include "separator.h"
#include "utils.h"
#include "TreeletEnumeratorAugmentedCycles.h"
using namespace std;
using namespace pandore;


void  MoleculeGraph::computeAugmentedCycleSpectrum(){
  //Nodes/Edges tables construction
  vector<string> nodes = this->_ccHypergraph_nodes_labels;
  vector<string> edges = this->_ccHypergraph_edges_labels;
  vector<std::map< std::pair<int,int>, unsigned short> > cycles(nodes.size());
  for(int i=0;i<_relevantCycles.size();i++)
    cycles[i] =  _relevantCycles[i]._angles;

  //Enumeration of treelets without hyperedges
  treelet_spectrum ** distrib = 
    TreeletEnumeratorAugmentedCycles::computeAugmentedTreeletSpectrum((*_ccHypergraph),nodes,edges,cycles);
  
  this->cc_augmented_spectrum = distrib;
  long cast_augmented_spectrum = (long)(this->cc_augmented_spectrum);
  _collection->SETVALUE("cc_augmented_spectrum", Long, (Long)(cast_augmented_spectrum));
}

treelet_spectrum **  MoleculeGraph::getAugmentedCycleSpectrum(){
    return cc_augmented_spectrum;
}

void MoleculeGraph::computeAugmentedCycles(){
  for(unsigned int i=0;i<_relevantCycles.size();i++){
    computeCycleSubs(i);
    _relevantCycles[i].computeAngles();
  }
#ifdef DEBUG

  for(unsigned int i=0;i<_relevantCycles.size();i++){
    cout << "Angles C_"<< i << " :" << endl;    
    Cycle c = _relevantCycles[i];
    map< std::pair<int,int>, unsigned short>::iterator it;
    for (it = c._angles.begin();it != c._angles.end();it++){
      cout << "Angle between " << it->first.first<<" and ";
      cout << it->first.second<<" : " << it->second << endl;
    }
  }
#endif
}


/*
** Returns a vector of size 2*nb_edges encoding the relative positioning of substituents around a cycle**
*/
vector<string> MoleculeGraph::computeCycleSubs(unsigned int node_index){

#ifdef DEBUG
  cout << "***** Cycle "<< node_index << ";"<< 
    MoleculeGraph::translateTreeletCycleCode(_relevantCycles[node_index]._code) <<  endl;
  cout << "Liste des substituants : " << endl;
#endif

  Cycle &c = _relevantCycles[node_index];

  vector<Substituent> subs(c._size *2);
  vector<string> subs_list;
  vector<string> subs_list_cool;
  
  vector<string> sub_vec_cool (c._size *2, ""); //vecteur des substituants
  vector<string> sub_vec (c._size *2, ""); //vecteur des substituants

  int * pos_sub = new int [_ccHypergraph_nodes_labels.size()]; /*position du sub correspondant 
								au noeud i dans subs_list*/
  memset(pos_sub,-1,_ccHypergraph_nodes_labels.size()*sizeof(int));
  
  // Parcours des voisins du cycle dans l'hypergraphe
  for(GEdge * e= (*_ccHypergraph)[node_index]->Neighbours();e!=NULL;e=e->Next()){
    stringstream  sub_label;
    sub_label << _ccHypergraph_edges_labels[e->Item()] <<"-";
    sub_label << _ccHypergraph_nodes_labels[e->Node()];
    subs_list.push_back(sub_label.str());
    
    stringstream sub_label_cool;
    sub_label_cool << translateTreeletCode(_ccHypergraph_edges_labels[e->Item()])<<"-";
    if(e->Node() >= _relevantCycles.size()) //Si le voisin *n'est pas* un cycle
      sub_label_cool <<	translateTreeletCode(_ccHypergraph_nodes_labels[e->Node()]);
    else //Si le voisin *est* un cycle
      sub_label_cool << translateTreeletCycleCode(_ccHypergraph_nodes_labels[e->Node()]);
    subs_list_cool.push_back(sub_label_cool.str());
    
    pos_sub[e->Node()] = subs_list.size()-1;
    
#ifdef DEBUG
    if(e->Node() >= _relevantCycles.size()) //Si le voisin *n'est pas* un cycle
      cout << _ccHypergraph_edges_labels[e->Item()]  << "-" << 
	translateTreeletCode(_ccHypergraph_nodes_labels[e->Node()]) << endl;
    else //Si le voisin *est* un cycle
      cout << translateTreeletCode(_ccHypergraph_edges_labels[e->Item()])  << "-" << 
	translateTreeletCycleCode(_ccHypergraph_nodes_labels[e->Node()]) << endl;
#endif
  }
  //subs_list contient l ensemble des substituants du cycle courant
  
#ifdef DEBUG
  // cout <<  "sub list : "<< endl;
  // for(int i =0;i<subs_list.size();i++)
  //   cout << i << " " << subs_list[i] << endl;
  cout <<  "sub list cool : "<< endl;
  for(int i =0;i<subs_list_cool.size();i++)
    cout << i << " " << subs_list_cool[i] << endl;
#endif
  
  //Identification des positions de chaque sub autour du cycle courant
  bool * edges_traversal = new bool[nbBonds];
  memcpy(edges_traversal, c._edges_print, sizeof(bool)*nbBonds); //copie des aretes du cycle c

  unsigned int start=-1;
  while(!c._nodes_print[++start]); //On trouve le 1er noeud de c
  unsigned int current = start;
  
  for(int i=0; i<c._size;i++){ //Parcours des c._size edges
    //Process node current
#ifdef DEBUG
    cout << "Iteration " << i << ", Noeud courant : "  << current << endl;
#endif
    for(GEdge * v= (*_graph)[current]->Neighbours();v!=NULL;v=v->Next()){
      if(! c._edges_print[v->Item()]){ //l'arete v n'appartient pas a c => substituant !
	Substituent s;
	s._position = i*2;
	for (int k=0;k<graph_to_hypergraph[v->Node()].size();k++){
	  s._hypergraph_node = graph_to_hypergraph[v->Node()][k];
	  if (pos_sub[graph_to_hypergraph[v->Node()][k]] != -1)
	    break;
	}
	if(pos_sub[s._hypergraph_node] == -1){
	}
	else{
	  
	string sub_label = subs_list[pos_sub[s._hypergraph_node]];
	s._code = sub_label;
	
	sub_vec[i*2] = sub_label;

	string sub_label_cool = subs_list_cool[pos_sub[s._hypergraph_node]];
	sub_vec_cool[i*2] = sub_label_cool;
		
#ifdef DEBUG
	s._code = sub_label_cool;
	cout << "Arete " << v->Item() << " : (" << current << "-"  << v->Node() << ")"   << " : " ;
	cout  << sub_label_cool << endl;
#endif
	subs[s._position] = s;
	}
      }
    }
    GEdge * e= (*_graph)[current]->Neighbours();
    while((e != NULL) && ! edges_traversal[e->Item()]){
      e=e->Next();
    }
    // Process edge e
    assert(e!=NULL); //Pb : un edge non trouvé.
    edges_traversal[e->Item()] = 0;
    current = e->Node();
  }
  //Raffinement des positions pour les cycles voisins
  vector<Substituent>::iterator it_sub;
  Substituent cur_sub = *(subs.begin());
  int positions =0;
  int nb_instances = 0;
  for(it_sub=subs.begin();it_sub!=subs.end();it_sub++){
    Substituent s = *it_sub;
    if(s._hypergraph_node >-1){
      if(s._hypergraph_node < _relevantCycles.size()){ //s correspond a un cycle
	int min_pos = s._position;
	int max_pos = s._position;
	vector<Substituent>::iterator loc_it_sub;
	for(loc_it_sub=subs.begin();loc_it_sub!=subs.end();loc_it_sub++){
	  if ((*loc_it_sub)._hypergraph_node == s._hypergraph_node){
	    min_pos = min((*loc_it_sub)._position,min_pos);
	    max_pos = max((*loc_it_sub)._position,max_pos);
	  }
	}
	int intersection_size = 0;
	Cycle adjacent = _relevantCycles[s._hypergraph_node];
	for(int i=0;i<nbAtoms();i++)
	  intersection_size += (c._nodes_print[i] && adjacent._nodes_print[i])?1:0;
	if ((max_pos-min_pos)%c._size == intersection_size)
	  s._position = (min_pos + (intersection_size/2))%(c._size*2);
	else
	  s._position = (max_pos + (intersection_size/2))%(c._size*2);
      }
      c._subs.insert(make_pair(s._hypergraph_node,s));
    }
  }

    
#ifdef DEBUG
  cout << endl;
  // cout << "Vecteur de sub :"<< endl;
  // for(int i =0;i<sub_vec.size();i++)
  //   cout << i << " " << sub_vec[i] << endl;
  cout << "Vecteur de sub cool :"<< endl;
  for(int i =0;i<sub_vec_cool.size();i++)
    cout << i << " " << sub_vec_cool[i] << endl;
  map<int, Substituent>::iterator it_map_sub;
  for ( it_map_sub=c._subs.begin(); it_map_sub!=c._subs.end(); ++it_map_sub){
    Substituent s = it_map_sub->second;
    cout << "Sub key "<< it_map_sub->first << " : position " << s._position << ", ";
    cout << s._code << endl;
  }
#endif
  return sub_vec;
}


void MoleculeGraph::Cycle::computeAngles(){
  map<int, Substituent >::iterator sub_it;
  for (sub_it = _subs.begin();sub_it != _subs.end();sub_it++){
    Substituent s1=sub_it->second;
    map<int, Substituent >::iterator sub_it2;
    for (sub_it2 = _subs.begin();sub_it2 != _subs.end();sub_it2++){
      Substituent s2=sub_it2->second;
      if(s1._hypergraph_node != s2._hypergraph_node){
	int n = _size*2;
	int i = s2._position - s1._position;
	unsigned short angle = (i % n + n) % n;
#ifdef DEBUG
	cout << "Angle between " << s1._hypergraph_node <<" and ";
	cout << s2._hypergraph_node <<" : " << angle << endl;
#endif
	_angles.insert(make_pair(make_pair(s1._hypergraph_node,s2._hypergraph_node),angle));
      }
    }
  }
}

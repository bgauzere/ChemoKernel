/*
 * @file graph.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Dec  5 2012
 * 
 * @todo 
 * @bug 
 *  
 * 
 * 
 *
 */

#include <iostream>
#include <map>
#include <sstream>
#include <algorithm>
#include "MoleculeGraph.h"
#include "string_utils.h"
#include "separator.h"
#include "utils.h"

using namespace std;
using namespace pandore;

void MoleculeGraph::computeContractedCycleGraphs(){
  this->computeRelevantCycles();
  if(_relevantCycles.size() == 0){ 
    _ccHypergraph = (Graph3d*)_graph->Clone();
    _contracted_ccHypergraph = NULL;
    //Copying labels
    for(int i=0;i<nbBonds;i++){
      stringstream s;
      s << _edges[i];
      _ccHypergraph_edges_labels.push_back(s.str());
    }
    for(int i=0;i<this->nbAtoms();i++){
      char c = _nodes[i]->GETVALUE("atom", Char);
      stringstream s;
      s << c;
      _ccHypergraph_nodes_labels.push_back(s.str());
    }
    /* No cycles to contract, we're done !*/
    return;
  }
  computeContractedCycleHypergraph();
  
  if(hyperedges.size() > 0)
    transformHypergraphToGraph();
  else
    _contracted_ccHypergraph = NULL;
}

void MoleculeGraph::computeContractedCycleHypergraph(){
  
  _ccHypergraph = (Graph3d*)_cyclicSystem->Clone();
#ifdef DEBUG
  cout << "Nb cyclic nodes : " << _ccHypergraph->Size() << endl;
#endif
  /*1st step : Include pure relevant cycle graph information*/
  //Copy nodes and edges labels
  for(vector<Cycle>::iterator it=_relevantCycles.begin();it<_relevantCycles.end();it++)
    _ccHypergraph_nodes_labels.push_back(it->_code);
  for(vector<string>::iterator it=_cyclicSystemBonds.begin();it<_cyclicSystemBonds.end();it++)
    _ccHypergraph_edges_labels.push_back(*it);
#ifdef DEBUG
  cout << "Nb cyclic edges : " <<  _ccHypergraph_edges_labels.size() << endl;
#endif
  
  /*2nd step : Include pure acyclic parts information*/
  /*Atoms*/
  // Encodes node projection from _graph to _ccHypergraph
  graph_to_hypergraph = new vector<int>[this->nbAtoms()];
  //Identification of nodes included into a relevant cycle 
  bool * atoms_in_cycles = new bool[this->nbAtoms()];
  memset(atoms_in_cycles,0,sizeof(bool)*this->nbAtoms());
  int nb_acyclic_atoms = 0;
  for(unsigned int j=0;j<_relevantCycles.size();j++)
    //vector<Cycle>::iterator it=_relevantCycles.begin();it<_relevantCycles.end();it++)
    for (int i=0;i<this->nbAtoms();i++)
      if(_relevantCycles[j]._nodes_print[i]){
	atoms_in_cycles[i] = true;
	graph_to_hypergraph[i].push_back(j);//g[i] fait partie du cycle _relevant_cycle[j],
	                                    //equivalent au noeud _ccHypergraph[j]
      }
  for (int i=0;i<this->nbAtoms();i++)
    nb_acyclic_atoms += (!atoms_in_cycles[i]);
  
#ifdef DEBUG
  for (int i=0;i<this->nbAtoms();i++)
    cout << atoms_in_cycles[i] << " ";
  cout << endl;
  cout << "Nb acyclic atoms : " << nb_acyclic_atoms << endl;
#endif

  //Adding acyclic atoms
  _ccHypergraph->Enlarge(_ccHypergraph->Size() + nb_acyclic_atoms);
  for (int i=0;i<this->nbAtoms();i++){
    if(!atoms_in_cycles[i]){
      char c = _nodes[i]->GETVALUE("atom", Char); //get atom label
      stringstream s;
      s << c;
      _ccHypergraph_nodes_labels.push_back(s.str());
      int index = _ccHypergraph_nodes_labels.size()-1;
      _ccHypergraph->Add(index,index);
      graph_to_hypergraph[i].push_back(index);//Normalement, graph_to_hypergraph[i] est vide
    }
  }
#ifdef DEBUG
  cout << "Nb nodes : " << _ccHypergraph->Size() << endl;
  cout << "Acyclic Edges : "  << endl;
#endif
  //Adding acyclic edges
  //NB : Convention : on connecte a et b si a < b
  int edge_index=_ccHypergraph_edges_labels.size();
  for(int i=0;i<_graph->Size();i++)
    if(!atoms_in_cycles[i])
      for(GEdge * e = (*_graph)[i]->Neighbours();e!=NULL;e=e->Next()){
	if (!atoms_in_cycles[e->Node()]){ //Edge between acyclic parts, we add it !
	  int n1 = graph_to_hypergraph[i][0];
	  int n2 = graph_to_hypergraph[e->Node()][0];
	  if(n1 < n2){
#ifdef DEBUG
	    cout << n1 << ", " <<  n2 << endl;
#endif
	    _ccHypergraph->Link(n1,n2,edge_index, 1.0f, false);
	    _ccHypergraph_edges_labels.push_back(string(1,_edges[e->Item()]));
	    edge_index ++;
	  }
	}
      }

  /*3rd step : Include other connections*/
  //Hyperedges storage:
  edges_to_hyperedges = new int[nbBonds]; //Correspondance egde de _graph et hyperedges
  memset(edges_to_hyperedges,-1,sizeof(int)*nbBonds);/*edges_to_hyperedges[i] is equals to -1 if
						       _edges[i] is not an hyperedge*/
  for(int v=0;v<_graph->Size();v++)
    if(atoms_in_cycles[v]){ // XXX: Factorisable avec precedent for
      bool multi_cycle_v = (graph_to_hypergraph[v].size() > 1); /*true if i belongs
								  to more than 1 cycle*/
#ifdef DEBUG
      //       cout << "Node " << v << " is ";
      //       if (multi_cycle_v)
      // 	cout << " multi cycle !"  << endl;
      //       else
      // 	cout << " not multi cycle !"  << endl;
#endif
	for(GEdge * e = (*_graph)[v]->Neighbours();e!=NULL;e=e->Next()){
	  int v_prime = e->Node();
	  //if(v < v_prime){ 
	    //Verif si e appartient a un edge d'un cycle de v dans le cas ou v appartient a un cycle
	    bool new_edge = true;
	    for(unsigned int i=0;i<graph_to_hypergraph[v].size();i++){
#ifdef DEBUG
	      cout << "Node " << v << " is in " << graph_to_hypergraph[v][i] << endl;
	      cout << "edge_print for cycle " << graph_to_hypergraph[v][i] << " :" ;
	      for(int j=0;j<this->nbBonds;j++)
		cout << _relevantCycles[graph_to_hypergraph[v][i]]._edges_print[j] << " "; 
	      cout << endl;
#endif
	      new_edge = new_edge && !_relevantCycles[graph_to_hypergraph[v][i]]._edges_print[e->Item()];
	      //False if e is included in one of same cycle as v.
	    }
	    if(new_edge){
	      bool multi_cycle_v_prime = (graph_to_hypergraph[v_prime].size() > 1);
	      if(multi_cycle_v || multi_cycle_v_prime){//Hyperarete
#ifdef DEBUG
		cout << "Hyperedge get ! : " << e->Item() << endl;
#endif
		if(edges_to_hyperedges[e->Item()] == -1)//Hyperedge not included, adding hyperedge
		  {
		    set<int> u(graph_to_hypergraph[v].begin(),graph_to_hypergraph[v].end());
		    set<int> v(graph_to_hypergraph[v_prime].begin(),graph_to_hypergraph[v_prime].end());
		    hyperedges.push_back(pair< set<int>, set<int> >(u,v));
		    edges_to_hyperedges[e->Item()] = hyperedges.size()-1;
		    hyperedges_to_edges.push_back(e->Item());
		  }
		else{//Update hyperedge
		  ;
#ifdef DEBUG
		  cout << "Need to update Hyperedge : Houston, we have a problem ...  " << e->Item() << endl;
#endif	       
		}
	      } else { //Classic edge
#ifdef DEBUG
		cout << "Classic edge get ! : " << e->Item()  << endl;
#endif
		//Testing if already added
		bool already_added = false;
		for(GEdge * e_2 = (*_ccHypergraph)[graph_to_hypergraph[v_prime][0]]->Neighbours();
		    e_2!=NULL;e_2=e_2->Next())
		  {
		    if(e_2->Node() == graph_to_hypergraph[v][0])
		      already_added = true;
		  }
		if(!already_added){
		  _ccHypergraph->Link(graph_to_hypergraph[v][0],
				      graph_to_hypergraph[v_prime][0],edge_index, 1.0f, false);
		  _ccHypergraph_edges_labels.push_back(string(1,_edges[e->Item()]));
		  edge_index ++;
		}
	      }
	    }
	    //}
	}
    }
#ifdef DEBUG
  cout << "Conctracted Cycle Hypergraph : "<< endl;
  cout << "Nb nodes : " <<  _ccHypergraph_nodes_labels.size() << endl;
  cout << "Nb edges : " <<  _ccHypergraph_edges_labels.size() << endl;
  cout << "Nb Hyperedges : " << hyperedges.size() << endl;
  
  cout << "Nodes : "<<endl;
  for(unsigned int i=0;i<_ccHypergraph_nodes_labels.size();i++)
    cout << "\t Node " << i << " : " << _ccHypergraph_nodes_labels[i] << endl;
  
  cout << "Edges : "<<endl;
  for(unsigned int i=0;i<_ccHypergraph_nodes_labels.size();i++)
    for(GEdge * e = (*_ccHypergraph)[i]->Neighbours();e!=NULL;e=e->Next())
      if((int)(i) < e->Node())
	cout << "\t Edge between " << i << " and  " << e->Node() << " : " << 
	  _ccHypergraph_edges_labels[e->Item()] << endl;
  
  cout << "Hyperedges : "<<endl;
  for(int e=0;e<nbBonds;e++)
    if(edges_to_hyperedges[e] > -1){
      cout << "Edge " << e << " is Hyperedge " << edges_to_hyperedges[e] << " between " << endl;
      cout << "\t u = ( " ;
      for(set<int>::iterator it = hyperedges[edges_to_hyperedges[e]].first.begin();
	  it != hyperedges[edges_to_hyperedges[e]].first.end();it++)
	cout << *it << ", ";
      cout << ")" << endl;
      
      cout << "\t v = ( " ;
      for(set<int>::iterator it = hyperedges[edges_to_hyperedges[e]].second.begin();
	  it != hyperedges[edges_to_hyperedges[e]].second.end();it++)
	cout << *it << ", ";
      cout << ")" << endl;
    }
  #endif

  delete [] atoms_in_cycles;
  //  delete [] graph_to_hypergraph;

  
}

int MoleculeGraph::contractSetofNodes(vector<int> vec_u){
  int node_u = *(vec_u.begin());
  if(vec_u.size() > 1){
    string label_u = computeContractedNodesCode(vec_u);
    for(vector<int>::iterator it =++(vec_u.begin());it!=vec_u.end();it++)
      {
	_contracted_ccHypergraph->Merge(node_u, *it);
	/*it est supprimé*/
	_contracted_ccHypergraph_nodes_labels[*it] = "";
      }
    _contracted_ccHypergraph_nodes_labels[node_u] = label_u;
  }
  return node_u;
}

void MoleculeGraph::transformHypergraphToGraph(){
  
  /*Hypergraph completed, contracting hyperedges*/
  _contracted_ccHypergraph = (Graph3d*)_ccHypergraph->Clone();
  _contracted_ccHypergraph_nodes_labels = _ccHypergraph_nodes_labels;
  _contracted_ccHypergraph_edges_labels = _ccHypergraph_edges_labels; 
  int current_nb_nodes = _contracted_ccHypergraph->Size();
  vector< vector<int> > classes;
  
  
  for(unsigned int e=0;e<hyperedges.size();e++){
    vector<int> v = vector<int>(hyperedges[e].first.begin(),hyperedges[e].first.end());
    sort(v.begin(),v.end());
    classes.push_back(v);
    v=vector<int>(hyperedges[e].second.begin(),hyperedges[e].second.end());
    sort(v.begin(),v.end());
    classes.push_back(v);
  }
  
  bool contraction=true;
  while(contraction){
    contraction = false;
    vector<int>::iterator it;
    for(unsigned int i=0;i<classes.size() && !contraction;i++)
      for(unsigned int j=i+1;j<classes.size() && !contraction;j++){
	vector<int>::iterator it;
	vector<int> inter = vector<int>(classes[i].size()+classes[j].size());
	if((it = set_intersection(classes[i].begin(),classes[i].end(),
				  classes[j].begin(),classes[j].end(),
				  inter.begin())) != inter.begin()){
	  //We contract !
	  vector<int> union_nodes = vector<int>(classes[i].size()+classes[j].size()-int(it - inter.begin()));
	  it = set_union(classes[i].begin(),classes[i].end(),
			 classes[j].begin(),classes[j].end(),union_nodes.begin());
	  classes[i] = union_nodes;
	  sort(classes[i].begin(),classes[i].end());
	  classes.erase(classes.begin()+j);
	  contraction = true;
	}	  
    }
  }
  //Classes ok
  vector<int> association_nodes_classes = vector<int>(current_nb_nodes);
  for(unsigned int i=0;i<classes.size();i++){
    //Contraction
    int index_classe = contractSetofNodes(classes[i]);
    for(unsigned int j=0;j<classes[i].size();j++){
      association_nodes_classes[classes[i][j]] = index_classe;
    }    
  }
  for(unsigned int e=0;e<hyperedges.size();e++){
    _contracted_ccHypergraph->Link(association_nodes_classes[*(hyperedges[e].first.begin())],
				   association_nodes_classes[*(hyperedges[e].second.begin())], 
				   _contracted_ccHypergraph_edges_labels.size(),1.0,false);
    stringstream label_edge;
    label_edge << _edges[hyperedges_to_edges[e]];
    _contracted_ccHypergraph_edges_labels.push_back(string(label_edge.str()));
  }
  

  
#ifdef DEBUG
  cout << "Hypergraphe contracté : "<< endl;
  cout << "Noeuds : "<< endl;
  for(int i=0;i<_contracted_ccHypergraph_nodes_labels.size();i++)
    if(_contracted_ccHypergraph_nodes_labels[i] != "")
      cout << i << " : " << _contracted_ccHypergraph_nodes_labels[i] << endl;
  
  cout << "Aretes : "<< endl;
  for(int i=0;i<_contracted_ccHypergraph_nodes_labels.size();i++)
    if(_contracted_ccHypergraph_nodes_labels[i] != "")
      for(GEdge * e = (*_contracted_ccHypergraph)[i]->Neighbours();e!=NULL;e=e->Next()){
	if(i < e->Node())
	  cout << i << ", " << e->Node() << " : " << _contracted_ccHypergraph_edges_labels[e->Item()] << endl;
      }
#endif
}


string MoleculeGraph::computeContractedNodesCode_rec(int current_node,bool * visited){
  visited[current_node] = true;
  bool all_visited = true;
  for(unsigned int i=0;i<_ccHypergraph_nodes_labels.size();i++)
    all_visited = all_visited && visited[i];
  
  if(all_visited){
    return _contracted_ccHypergraph_nodes_labels[current_node];
  }
  else /*On calcule tous les codes possibles en partant de la*/
    {
      stringstream s;
      s << _ccHypergraph_nodes_labels[current_node];
      s << CONTRACTED_NODE_SEP;
      vector<string> codes;
      for(GEdge * e = (*_ccHypergraph)[current_node]->Neighbours();e!=NULL;e=e->Next())
	{
	  if(!visited[e->Node()]){//if not visited, v->Node \in set_of_nodes
	    stringstream s_prime;
	    s_prime << s.str();
      	    s_prime << _ccHypergraph_edges_labels[e->Item()];
      	    s_prime << CONTRACTED_NODE_SEP;
	    bool * visited_cur = new bool[_ccHypergraph_nodes_labels.size()];
	    memcpy(visited_cur,visited,sizeof(bool)*_ccHypergraph_nodes_labels.size());
	    string new_code = computeContractedNodesCode_rec(e->Node(),visited_cur);
	    if(new_code != ""){	      //Candidat iff all nodes have been visited.
	      s_prime << new_code;
	      codes.push_back(s_prime.str());
	    }
	    delete [] visited_cur;
	  }
	}
      
      if(codes.size() > 0){
	sort(codes.begin(),codes.end());
	return codes[0];
      }else
	return "";
    }
}

string MoleculeGraph::computeContractedNodesCode(vector<int> set_of_nodes){
  //set_of_nodes is a subset of cc_Hypergraph
  //We have to find the spanning tree over set_of_nodes having minimal code according to depth traversal

  vector<string> codes;
  bool * visited = new bool[_ccHypergraph_nodes_labels.size()];
  memset(visited,true,sizeof(bool)*_ccHypergraph_nodes_labels.size());
  for(unsigned int i=0;i<set_of_nodes.size();i++)
    visited[set_of_nodes[i]] = false;
  
  for(unsigned int i =0; i<set_of_nodes.size();i++){
    bool * visited_cur = new bool[_ccHypergraph_nodes_labels.size()];
    memcpy(visited_cur,visited,sizeof(bool)*_ccHypergraph_nodes_labels.size());
    string code = computeContractedNodesCode_rec(set_of_nodes[i],visited_cur);
    if(code != "")
      codes.push_back(code);
    delete [] visited_cur;
  }
  delete [] visited;
  assert(codes.size() > 0);
  sort(codes.begin(),codes.end());
  return codes[0];
}

void MoleculeGraph::computeContractedCycleSpectrum(){
  //Nodes/Edges tables construction
  vector<string> nodes = this->_ccHypergraph_nodes_labels;
  vector<string> edges = this->_ccHypergraph_edges_labels;

  //Enumeration of treelets without hyperedges
  treelet_spectrum ** distrib = TreeletEnumerator::computeTreeletSpectrum((*_ccHypergraph),nodes,edges);
  assert(distrib != NULL);
  //Enumeration of treelets with hyperedges
  if(_contracted_ccHypergraph != NULL){
    nodes = this->_contracted_ccHypergraph_nodes_labels;
    edges = this->_contracted_ccHypergraph_edges_labels;
    treelet_spectrum ** distrib_2 = TreeletEnumerator::computeTreeletSpectrum((*_contracted_ccHypergraph),nodes,edges);
    //Union of two treelets distributions
    for(int i=0;i<SIZE_SPECTRUM;i++){
      for(map<string,double>::iterator it=distrib_2[i]->begin();it!=distrib_2[i]->end();it++)
	if(it->first.find(CONTRACTED_NODE_SEP_CHAR) > std::string::npos) //We only keep treelets involving ContractedNodes 
	  distrib[i]->insert(pair<string, double>(it->first,it->second));
    }
  }
  this->cc_hypergraph_spectrum = distrib;
  long cast_contracted_cycle_spectrum = (long)(this->cc_hypergraph_spectrum);
  _collection->SETVALUE("cc_hypergraph_spectrum", Long, (Long)(cast_contracted_cycle_spectrum));
}

treelet_spectrum ** MoleculeGraph::getContractedCycleSpectrum(){  /*The spectrum must be computed*/ 
  return cc_hypergraph_spectrum;
}




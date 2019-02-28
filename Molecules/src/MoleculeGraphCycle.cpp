/*
 * @file MoleculeGraphCycle.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Oct 18 2011
 * 
 * @todo
 * @bug 
 *  
 * Methods for cycle enumeration of a molecular graph
 * 
 *
 */
#include "MoleculeGraph.h"

#include <iostream>
#include <stack>
#include <algorithm>
#include <limits.h>
#include <assert.h>

#include "utils.h"
#include "separator.h" 

using namespace std;
using namespace pandore;
using namespace cimg_library;

int * numbering; /*atoms numbering of Tarjan algorithm*/
int * lowpt; /*lowpt of Tarjan algorithm*/
int i;
stack< pair<int,int> > edge_stack;

void MoleculeGraph::computeTwoconnexeComponents()
{
  /*Initializations*/
  i=1;
  numbering = new int[nbAtoms()];
  memset(numbering, 0,nbAtoms()*sizeof(int));
  lowpt = new int[nbAtoms()];
  memset(lowpt, 0,nbAtoms()*sizeof(int));
  edge_stack = stack< pair<int,int> >();
  _twoconnexeComponents = vector< vector< pair <int, int> > >();
  for(int w=0; w<nbAtoms(); w++)
    {
      if(numbering[w] == 0) //w isn't numbered
  	biconnect(w,-1);
    }
}

void MoleculeGraph::biconnect(int v, int u)
{
  numbering[v] = i;
  i++;
  lowpt[v] = numbering[v];
  for(GEdge * w = (*_graph)[v]->Neighbours();w!=NULL;w=w->Next())
    {
      if(numbering[w->Node()] == 0) //w isn't numbered
	{
	  edge_stack.push( pair<int,int>(v,w->Node()));
	  biconnect(w->Node(),v);
	  lowpt[v] = min(lowpt[v], lowpt[w->Node()]);
	  if (lowpt[w->Node()] >= numbering[v])
	    {
	      /*New Biconnected component*/
	      vector< pair<int, int > > current_component;
	      stack< pair<int,int> > edge_stack_new = stack< pair<int,int> >();
	      while(! edge_stack.empty())
		{
		  pair<int,int> cur_edge = edge_stack.top();
		  edge_stack.pop();
		  if(numbering[cur_edge.first] > numbering[v])
		    {
		      current_component.push_back(cur_edge);
		    }
		  else if((cur_edge.first == v) && (cur_edge.second == w->Node()))
		    {
		      current_component.push_back(cur_edge);
		    }
		  else
		    {
		      edge_stack_new.push(cur_edge);
		    }
		}
	      _twoconnexeComponents.push_back(current_component);
	      edge_stack = edge_stack_new;
	    }
	}
      else if((numbering[w->Node()] < numbering[v]) && (w->Node() != u))
	{
	  edge_stack.push(pair<int,int>(v,w->Node()));
	  lowpt[v] = min(lowpt[v], numbering[w->Node()]);
	}
    }
}



/*Vismara 2003*/
int MoleculeGraph::getAncestor(int node, int root, vector<int> * parent)
{
  int current_node = node;
  while(! parent[current_node].empty())
    {
      if(parent[current_node][0] ==root)
	return current_node;
      current_node = parent[current_node][0];
    }
  
  //Error
  return -1;
}



void MoleculeGraph::computeVismaraCycles()
{
  // _cyclesCodes = vector<string>();
  vector<MoleculeGraph::Cycle> cycles;
  for(int r=0; r<nbAtoms(); r++)
    {
      int * levels = new int[nbAtoms()]; 
      for(int i=0;i<nbAtoms();i++)
	levels[i] = INT_MAX;
      levels[r] = 0; 
      vector<int> * parent = new vector<int>[nbAtoms()]; parent[r].push_back(r);

      vector<int> * parent_edges = new vector<int>[nbAtoms()]; //Pour la construction du code
      
      vector<int> * equals = new vector<int>[nbAtoms()];
      vector<int> visited;
      vector<int> * sphere = new vector<int>();sphere->push_back(r);
      vector<int> Xr; 
      for(int i =r; i <nbAtoms();i++)
	Xr.push_back(i);

      vector<int> v(nbAtoms());
      vector<int>::iterator it; 
      // it = set_intersection(Xr.begin(),Xr.end(), sphere->begin(), sphere->end(), v.begin());
      // cout << (it - v.begin()) << endl;
      sort(sphere->begin(), sphere->end());
      
      while (((it = set_intersection(Xr.begin(),Xr.end(), 
				     sphere->begin(), sphere->end(), v.begin())) - v.begin()) != 0)
      	{
	  vector<int> * copie_sphere = sphere;
      	  sphere = new vector<int>();
	  // for(int index_x = 0; index_x < copie_sphere->size(); index_x ++)
	  //   cout << copie_sphere->at(index_x) << ", ";
	  // cout << endl;
	  for(unsigned int index_x = 0; index_x < copie_sphere->size(); index_x ++)
	    {
	      int x = copie_sphere->at(index_x);
	      for(GEdge * e = (*_graph)[x]->Neighbours();e!=NULL;e=e->Next())
		{
		  int y = e->Node();
		  if(levels[x] < levels[y])
		    {
		      if ((y >= r) && (! parent[x].empty()))
			{
			  parent[y].push_back(x);
			  parent_edges[y].push_back(e->Item());
			}
		      levels[y] = levels[x] + 1;
		      bool find=false;
		      for(unsigned int i=0;i<sphere->size();i++)
			if (sphere->at(i) ==  y)
			  find = true;
		      if(!find)
			sphere->push_back(y);
		    }
		  else if ((!parent[y].empty()) && (!parent[x].empty()))
		    {
		      equals[y].push_back(x);
		    }
		  
		}
	      //On découvre les cycles
	      if(parent[x].size() >= 2)/*Cycle pair*/
		{
		  for(unsigned int index_p=0;index_p<parent[x].size();index_p++)
		    for(unsigned int index_q=index_p;index_q<parent[x].size();index_q++)
		      {
			int p = parent[x][index_p];
			int q = parent[x][index_q];
			
			int edge_p = parent_edges[x][index_p];
			int edge_q = parent_edges[x][index_q];

			int anc_p = getAncestor(p,r,parent);
			int anc_q = getAncestor(q,r,parent);
			
			if ((p != q) && (anc_p != -1) && 
			    (anc_q != -1) && (anc_p != anc_q))
			  {
			    vector<int> current_cycle;
			    // current_cycle.insert(current_cycle.begin(),pair<int,int>(x,p));
			    // current_cycle.push_back(pair<int,int>(x,q));
			    current_cycle.insert(current_cycle.begin(),edge_p);
			    current_cycle.push_back(edge_q);
			    int current_node = p;
			    while(current_node!= r)
			      {
			    	current_cycle.insert(current_cycle.begin(),
			    			     parent_edges[current_node][0]);
			    	current_node = parent[current_node][0];
			      }
			    current_node = q;
			    while(current_node!= r)
			      {
			    	current_cycle.push_back(parent_edges[current_node][0]);
			    	current_node = parent[current_node][0];
			      }
			    //Print Cycle
			    // cout << "cycle pair:";
			    // for (unsigned int i=0;i<current_cycle.size();i++)
			    //   cout << "(" << _edgesIncidentNodes[current_cycle[i]].first << "," 
			    // 	   <<  _edgesIncidentNodes[current_cycle[i]].second << "),";
			    // cout << endl;
			    recordCycle(current_cycle);
			  }
		      }
		}
	      
	      /*Cycle impair*/
	      sort(equals[x].begin(), equals[x].end());
	      sort(visited.begin(), visited.end());
	      vector<int> joints(nbAtoms());
	      vector<int>::iterator it_joints; 
	      if(((it_joints = set_intersection(equals[x].begin(),equals[x].end(), 
						visited.begin(), visited.end(), 
						joints.begin())) - joints.begin()) != 0)
		{
		  vector<int>::iterator it_e;
		  for(it_e=joints.begin();it_e!=it_joints;it_e++)
		    {		      
		      if(getAncestor(x,r,parent) != (getAncestor(*it_e,r,parent)))
			{
			  // cout << "cycle impair: " ;
			  vector<int> current_cycle;
			  int current_node = x;
			  for(GEdge * e = (*_graph)[x]->Neighbours();e!=NULL;e=e->Next())
			    {
			      if(e->Node() == *it_e) 
				current_cycle.push_back(e->Item());
			    }
			  while(current_node!= r)
			    {
			      current_cycle.push_back(parent_edges[current_node][0]);
			      current_node = parent[current_node][0];
			    }
			  current_node = *it_e;
			  while(current_node!= r)
			    {
			      current_cycle.insert(current_cycle.begin(),parent_edges[current_node][0]);
			      current_node = parent[current_node][0];
			    }
			  //Print Cycle
			  // for (unsigned int i=0;i<current_cycle.size();i++)
			  //   cout << "(" << _edgesIncidentNodes[current_cycle[i]].first << "," 
			  // 	 <<  _edgesIncidentNodes[current_cycle[i]].second << "),";
			  // cout << endl;
			  recordCycle(current_cycle);
			}
		    }
		}
	      
	      visited.push_back(x);
	      // cout << "visited: ";
	      // for(int i =0;i<visited.size();i++)
	      // 	cout << visited[i] << ",";
	      // cout << endl;
	    }
	  delete copie_sphere;
      	}
      delete [] parent;
      delete [] parent_edges;
      delete [] equals;
      delete sphere;
      delete [] levels;
    }
}


void MoleculeGraph::computeRelevantCycles()
{
  //Computing CI
  computeVismaraCycles();
  //CI Sorting
  // cout << "Before sort:" << endl;
  // for(int i=0;i<_cycles.size();i++)
  //   cout << _cycles[i]._code << endl;
  if(_cycles.size() == 0)
    return;
  sort(_cycles.begin(),_cycles.end(), cmp_cycle);
    
  // cout << "After sort:" << endl;
  // for(int i=0;i<_cycles.size();i++)
  //   cout << _cycles[i]._code << endl;
  
  ///*Prints cycles*/
  // for(int i=0;i<_cycles.size();i++)
  //   {
  //     for(int e=0;e<nbBonds;e++)
  // 	cout << _cycles[i]._edges_print[e] << " ";
  //     cout << endl;
  //   }
  //Calcul minimal basis
  CImg<bool> mat;
  
  // //Rajout d'un cycle pour test
  // Cycle test;
  // test._edges_print = new bool[nbBonds];
  // memset(test._edges_print,true,sizeof(bool)*nbBonds);
  // test._edges_print[nbBonds-1] = false;
  // test._code = string("cycle test");
  // _cycles.push_back(test);
  
  //Print selected_cycles
  // cout << "All Cycles:"<< endl;
  // for(unsigned int i=0;i<_cycles.size();i++)
  //   cout << _cycles[i]._code << endl;

  int current_size =_cycles[0]._size;
  int nb_cycles =0;
  for(unsigned int i = 0; i< _cycles.size(); i++)
    {
      CImg<bool>tmp_mat = mat;
      tmp_mat.resize(nbBonds,nb_cycles+1,1,1,0);
      
      for(int e = 0; e < nbBonds; e ++)
	tmp_mat(e,nb_cycles) = _cycles[i]._edges_print[e];
      //cout << "h:" << tmp_mat.height() <<",rang = " << utils::booleanRank(tmp_mat) << endl;
      // cout << "Mat:" << endl;
      // for(int y=0;y<tmp_mat.height();y++){
      // 	for(int x=0;x<tmp_mat.width();x++)
      // 	  cout << tmp_mat(x,y) << " " ;
      // 	cout << endl;
      // }


      if(utils::booleanRank(tmp_mat) == nb_cycles+1)/*C is relevant*/
	{
	  // cout << "Ajout" << endl;
	  mat = tmp_mat;
	  // cout << "Mat modif :" << endl;
	  // for(int y=0;y<mat.height();y++){
	  //   for(int x=0;x<mat.width();x++)
	  //     cout << mat(x,y) << " " ;
	  //   cout << endl;
	  // }

	  _relevantCycles.push_back(_cycles[i]);
	  nb_cycles++;
	}
      if((nb_cycles >= getCyclomaticNumber()) &&
	 (i+1 < _cycles.size()) &&
	 (_cycles[i+1]._size > current_size))
	break;
    }
  
  // cout << getCyclomaticNumber() << endl;
  // cout << nb_cycles << endl;
  
  //Print selected_cycles
  // cout << "Relevant Cycles:"<< endl;
  // for(unsigned int i=0;i<_relevantCycles.size();i++)
  //   cout << _relevantCycles[i]._code << endl;
  
  computeRelevantCycleGraph();
}

/*
** Returns canonical code for cycle**
*/

string MoleculeGraph::computeCycleCode(vector<int> c)
{
  string code;
  //Looking for the smallest 
  unsigned char min_label = 255;
  vector<int> min_positions;
  //Assume that bonds are continuous
  //Reorganize bonds
  for(unsigned int i=0;i<c.size();i++)
    {
      //Find the common node
      if ((_edgesIncidentNodes[c[i]].first == _edgesIncidentNodes[c[(i+1)%c.size()]].first) ||
  	  (_edgesIncidentNodes[c[i]].first == _edgesIncidentNodes[c[(i+1)%c.size()]].second))
  	{
	  //Reorganisation
  	  int temp = _edgesIncidentNodes[c[i]].first;
	  _edgesIncidentNodes[c[i]].first = _edgesIncidentNodes[c[i]].second;
  	  _edgesIncidentNodes[c[i]].second = temp;
  	}
      const unsigned char label = _nodes[_edgesIncidentNodes[c[i]].first]->GETVALUE("atom", Char);
      code.append(1,label);
      if(label <= min_label)
	{
	  if(label < min_label)
	    {
	      min_label = label;
	      min_positions.clear();
	    }
	  min_positions.push_back(i*2);//Code : NodeBondNodeBond...
	}
      // code.append(CYCLE_LABEL_SEP);
      char bond;
      bond = _edges[c[i]];
      // code.append(1,bond);
      code.append(1,'1');
      
      // code.append(CYCLE_LABEL_SEP);
    }
  //Minimize code
  vector<string> permutations;
  for(unsigned int i=0;i<min_positions.size();i++)
    {
      string permutation(code.begin() + min_positions[i], code.end());
      permutation.append(code.begin(), code.begin() + min_positions[i]);
      //permutation computed, adding the first node to the end for reverse code comparaison
      permutation.append(1,permutation[0]);
      permutations.push_back(permutation);
      permutations.push_back(utils::reverse(permutation));
    } 
  //Sort them
  sort(permutations.begin(),permutations.end());
  return permutations[0];
}


/*
** Returns canonical code of a link between two cycles**
*/
string MoleculeGraph::computeCycleLinkCode(vector<int> c)
{
  string code;
  //Assume that bonds are continuous
  //Reorganize bonds
  for(unsigned int i=0;i<c.size();i++)
    {
      //Find the common node
      if ((i != c.size()-1) && //Next doesn't exist !
	  ((_edgesIncidentNodes[c[i]].first == _edgesIncidentNodes[c[i+1]].first) ||
	   (_edgesIncidentNodes[c[i]].first == _edgesIncidentNodes[c[i+1]].second)))
	{ //First node belongs to c[i] and c[i+1]
	  //Swapping with second node
	  //Reorganisation
  	  int temp = _edgesIncidentNodes[c[i]].first;
	  _edgesIncidentNodes[c[i]].first = _edgesIncidentNodes[c[i]].second;
  	  _edgesIncidentNodes[c[i]].second = temp;
  	}
      const char label = _nodes[_edgesIncidentNodes[c[i]].first]->GETVALUE("atom", Char);
      code.append(1,label);
      string bond;
      bond = _edges[c[i]];
      code.append(bond);
    }
  //Adding end of sequence
  const char label = _nodes[_edgesIncidentNodes[c[c.size()-1]].second]->GETVALUE("atom", Char);
  code.append(1,label);
  //minimize code
  code = (code>utils::reverse(code))?utils::reverse(code):code;
  return code;
}

void MoleculeGraph::recordCycle(vector<int> c)
{
  Cycle current_cycle;
  current_cycle._edges_print = new bool[nbBonds];
  current_cycle._nodes_print = new bool[nbAtoms()];

  memset(current_cycle._edges_print, 0, sizeof(bool)*nbBonds);
  memset(current_cycle._nodes_print, 0, sizeof(bool)*nbAtoms());
  
  for(unsigned int i=0;i<c.size();i++){
    current_cycle._edges_print[c[i]] = true;
    current_cycle._nodes_print[_edgesIncidentNodes[c[i]].first] = true;
    current_cycle._nodes_print[_edgesIncidentNodes[c[i]].second] = true;
  }
  
  current_cycle._code = computeCycleCode(c);
  current_cycle._size = c.size();
  _cycles.push_back(current_cycle);
}

void MoleculeGraph::computeRelevantCyclesGraphsBonds(){
  int nb_edges = 0;
  
  for(unsigned int i =0; i<_relevantCycles.size();i++)
    for(unsigned int j=i+1; j<_relevantCycles.size();j++)
      {
	vector<int> common_nodes;
	vector<int> common_edges;
	int nb_comp = 0;
  	//Test connection between cycles i and j
  	//Cv(i) \inter Cv(j)
  	// cout << "Intersection des noeuds: " << endl;
  	for(int n=0;n<nbAtoms();n++)
  	  {
  	    if (_relevantCycles[i]._nodes_print[n] && _relevantCycles[j]._nodes_print[n])
  	      {
		common_nodes.push_back(n);
  		// cout << n  << " : (" << (int)_nodes[n]->GETVALUE("atom", Char) << ")" << endl;
  	      }
  	  }
	if(common_nodes.size()== 1)
	  {
	    _cyclicSystemBonds.push_back(string(1,_nodes[common_nodes[0]]->GETVALUE("atom", Char)));
	    _cyclicSystem->Link((Long)(i),(Long)(j),(Long)(nb_edges++));
	  }
	else if(common_nodes.size() > 1)
	  {
	    //Ce(i) \inter Ce(j)
	    // cout << "Intersection des aretes: " << endl;
	    for(int n=0;n<nbBonds;n++)
	      {
		if (_relevantCycles[i]._edges_print[n] && _relevantCycles[j]._edges_print[n])
		  {
		    common_edges.push_back(n);
		    // cout << n << " : (" << _edges[n] << ")" << endl;
		  }
	      }
	    if(common_edges.size() >= 1) //XXX: Combinable avec etape précédente ?
	      {
		//Continuity detection:Seed propagation
		// mu = |N| - |E| + nb_comp, mu = 0 =>
		nb_comp = common_nodes.size() - common_edges.size();
		bool * visited = new bool[nbBonds];
		memset(visited,0,sizeof(bool)*nbBonds);
		for(unsigned int n=0;n<common_edges.size();n++)
		  if(!visited[common_edges[n]] )
		    {
		      vector<int> connexe_edge = getConnexeEdgeSequence(common_edges[n], common_edges, visited);
		      //Construction Code
		      _cyclicSystemBonds.push_back(computeCycleLinkCode(connexe_edge));
		      //_cyclicSystemBonds.push_back("");//XXX: Modif pour test
		      _cyclicSystem->Link((Long)(i),(Long)(j),(Long)(nb_edges++));
		    }
		delete [] visited;
	      }
	  }
      }
}


vector<int> MoleculeGraph::getConnexeEdgeSequence(int seed, vector<int> edges, bool * visited){
  vector<int> edge_sequence;
  edge_sequence.push_back(seed);
  visited[seed] = true;
  //Left visit
  int left_node = _edgesIncidentNodes[seed].first;
  for(GEdge * e = (*_graph)[left_node]->Neighbours();e!=NULL;e=e->Next())
    {
      bool is_e_common_edge = (find(edges.begin(),edges.end(), e->Item()) != edges.end());
      if(is_e_common_edge && !visited[e->Item()])
	{
	  vector<int> left_seq = getConnexeEdgeSequence(e->Item(), edges, visited);
	  edge_sequence.insert(edge_sequence.begin(), left_seq.begin(), left_seq.end());
	}      
    }

  
  //Right visit
  int right_node = _edgesIncidentNodes[seed].second;
  for(GEdge * e = (*_graph)[right_node]->Neighbours();e!=NULL;e=e->Next())
    {
      bool is_e_common_edge = (find(edges.begin(),edges.end(), e->Item()) != edges.end());
      if(is_e_common_edge && !visited[e->Item()])
	{
	  vector<int> right_seq = getConnexeEdgeSequence(e->Item(), edges, visited);
	  edge_sequence.insert(edge_sequence.end(), right_seq.begin(), right_seq.end());
	}      
    }
  return edge_sequence ;
}


void MoleculeGraph::computeRelevantCycleGraph(){
  
  _cyclicSystem = new Graph3d(false);
  _cyclicSystem->New(_relevantCycles.size(), 0, 0, 0);
  for(unsigned int i =0; i<_relevantCycles.size();i++)
    _cyclicSystem->Add(i,i);
  computeRelevantCyclesGraphsBonds();
  //Print Cyclic System
  // cout << "Cyclic System:"<< endl;
  // cout << " \tNodes:"<< endl;
  // for(int i=0;i<_cyclicSystem->Size();i++)
  //   cout << i << " : " << translateTreeletCode(_relevantCycles[(*_cyclicSystem)[i]->Item()]._code) << endl;
  // cout << " \tEdges:"<< endl;  
  // for(int i=0;i<_cyclicSystem->Size();i++)
  //   for(GEdge * e = (*_cyclicSystem)[i]->Neighbours();e!=NULL;e=e->Next())
  //     cout << i << " is connected to  " << e->Node() << " : " <<  translateTreeletCode(_cyclicSystemBonds[e->Item()]) << endl;
}

void MoleculeGraph::computeCycleTreeletSpectrum(){
  //Nodes/Edges tables construction
  if(this->_relevantCycles.size() == 0){    //Molecule acyclique
    this->cycle_spectrum = new treelet_spectrum*[SIZE_SPECTRUM];
    for(int i=0;i<SIZE_SPECTRUM;i++)
      this->cycle_spectrum[i] = new treelet_spectrum();
  }else{
    vector<string> nodes(this->_relevantCycles.size());
    for(unsigned int i=0;i<this->_relevantCycles.size();i++){
      nodes[i] = string(_relevantCycles[i]._code);
    }
    
    vector<string> edges(_cyclicSystemBonds.size());
    for(unsigned int i=0;i<_cyclicSystemBonds.size();i++){
      edges[i] = string(_cyclicSystemBonds[i]);
    }
    treelet_spectrum ** distrib = TreeletEnumerator::computeTreeletSpectrum((*_cyclicSystem),nodes,edges);
    assert(distrib != NULL);
    this->cycle_spectrum = distrib;
  }
  long cast_cycle_spectrum = (long)(this->cycle_spectrum);
  _collection->SETVALUE("cycle_spectrum", Long, (Long)(cast_cycle_spectrum));
}
treelet_spectrum ** MoleculeGraph::getCycleTreeletSpectrum(){  /*The spectrum must be computed*/ 
  return cycle_spectrum;
}

void MoleculeGraph::computeSimpleCycles(int k){
  vector< vector<Cycle> > _simpleCycles_iter(k);
  
  for (unsigned int i=0;i<_relevantCycles.size();i++){
    _simpleCycles_iter[0].push_back(_relevantCycles[i]);
    _simpleCycles.push_back(_relevantCycles[i]);
  }
  // if(_simpleCycles.size() ==0)
  //   return;
  for(int n=1;n<k;n++){
    vector<Cycle> cycles_current_iter;
    for (unsigned int i=0;i<_simpleCycles_iter[0].size();i++)
      for (unsigned int j=0;j<_simpleCycles_iter[n-1].size();j++){
	// cout << i << "," << j << "," << n << endl;
	Cycle tmp;
	tmp._edges_print = new bool[nbBonds];
	memset(tmp._edges_print, 0,this->nbBonds*sizeof(bool));
	tmp._nodes_print = new bool[nbAtoms()];
	memset(tmp._nodes_print, 0,this->nbAtoms()*sizeof(bool));
	vector<int> edges;
	for(int e=0;e<this->nbBonds;e++){
	  tmp._edges_print[e] = _simpleCycles_iter[0][i]._edges_print[e] ^ _simpleCycles_iter[n-1][j]._edges_print[e];
	  if(tmp._edges_print[e])
	    edges.push_back(e);
	}
	
	for(unsigned int e=0;e<edges.size();e++){
	  tmp._nodes_print[_edgesIncidentNodes[edges[e]].first] = true;
	  tmp._nodes_print[_edgesIncidentNodes[edges[e]].second] = true;
	}
	tmp._size = edges.size();
	if(tmp._size != 0){
	  //Test simple cycle
	  int start = _edgesIncidentNodes[edges[0]].first;
	  int current = 0;
	  int next = _edgesIncidentNodes[edges[0]].second;
	  int nb_edges = 1;
	  bool * visited_edges = new bool[edges.size()];
	  memset(visited_edges,0,edges.size()*sizeof(bool));
	  while(next != start){
	    for(int e=0;e<(int)edges.size();e++){
	      if((e!=current) && (((_edgesIncidentNodes[edges[e]].first == next) 
				   && (!visited_edges[e])) || 
				  ((_edgesIncidentNodes[edges[e]].second == next) 
				   && (!visited_edges[e])))){
		current=e;
		visited_edges[e] = true;
		next = (_edgesIncidentNodes[edges[e]].first == next)?
		  _edgesIncidentNodes[edges[e]].second:
		  _edgesIncidentNodes[edges[e]].first;

		e = edges.size();
		nb_edges ++;
	      }
	    }
	  }
	  if(nb_edges != tmp._size){
	    // cout << "Not a simple cycle : Only " << nb_edges << "/" << tmp._size << endl;
	  }else{
	    //Test if already found
	    bool found = false;
	    for(unsigned int c=0;c<_simpleCycles.size();c++){
	      bool same = true;
	      for(int e=0;e<this->nbBonds;e++)
		same = same & (_simpleCycles[c]._edges_print[e] == tmp._edges_print[e]);
	      found = found | same;
	    } 
	    if(!found){
	      // cout << "Cycle level " << n << endl;
	      tmp._code = computeCycleCode(edges);
	      cycles_current_iter.push_back(tmp);
	      _simpleCycles.push_back(tmp);
	    }
	  }
	}
      }
    for(unsigned int i=0;i<cycles_current_iter.size();i++)
      _simpleCycles_iter[n].push_back(cycles_current_iter[i]);
  }

  for(unsigned int i=0;i<_simpleCycles.size();i++){
      _simpleCycleCodes.push_back(_simpleCycles[i]._code);
  }

  long cast_simple_cycle = (long)(&(this->_simpleCycleCodes));
  _collection->SETVALUE("simple_cycle_codes", Long, (Long)(cast_simple_cycle));
}

vector<string >  MoleculeGraph::getSimpleCycleCodes(){
  return _simpleCycleCodes;
}

bool MoleculeGraph::isCyclicNode(int graph_node){
  for(int i = 0; i<_cycles.size();i++)
    if (_cycles[i]._nodes_print[graph_node])
      return true;
  return false;

}
bool MoleculeGraph::isCyclicEdge(int graph_edge){
}

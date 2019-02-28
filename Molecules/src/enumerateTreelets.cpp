/*
 * @file enumerateTreelets.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue May 15 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include "pandore.h"
#include "TreeletEnumerator.h"
using namespace std;
using namespace pandore;

vector<char*> split (const char* chaine, const char* sep)
{
  vector<char*> v;
	
  char* s = strtok ((char*)chaine, (char*)sep);
	
  while (s != NULL)
    {
      v.push_back (s);
      s = strtok (NULL, (char*)sep);
    }
	
  return v;
}


Graph3d* createGraph(char * filename, vector<string>& nodes,vector<string>& edges){
  int nbBonds;
  int nbAtoms;
	
  ifstream file(filename,ios::in);
  vector<char*> v;
  Graph3d * graph;
  
  if (file) 
    {		
      char * s = new char[255];
		
      file.getline(s, 255); // The first line is useless
      file.getline(s, 255); // Second line = NumberOfAtoms NumberOfBonds
      
      
      v = split(s, " ");
      
      nbAtoms = atoi(v[0]);
      nbBonds = atoi(v[1]);
		
      graph = new Graph3d(false);
      nodes = vector<string>(nbAtoms); 
      edges = vector<string>(nbBonds); 
      pair<int,int> * _edgesIncidentNodes = new pair<int,int>[nbBonds];

      // Graph initialization
      graph->New(nbAtoms, 0, 0, 0);
		
      // Creation of the nodes
		
      for(int i=0; i<nbAtoms; i++)
	{
	  file.getline(s,255); // s = x y z Label
	  v = split(s," ");
	  
	  string label = string(&(v[3][0]));
	  nodes[i] = label;
	  graph->Add(i,i); 
	}
		
      // Creation of the edges
      for (int i=0; i<nbBonds; i++)
	{
	  file.getline(s,255); // s = Node1 Node2 Label Label
	  v = split(s," ");
	  edges[i] = v[2];
	  
	  _edgesIncidentNodes[i].first = min(atoi(v[0])-1, atoi(v[1])-1);
	  _edgesIncidentNodes[i].second = max(atoi(v[0])-1, atoi(v[1])-1);
	  graph->Link(_edgesIncidentNodes[i].first, _edgesIncidentNodes[i].second,i,1.0f,false);
	}
		
      file.close();
      return graph;
    }
  else
    {
      cerr << "Error : Can't open file : " << filename << endl;
      graph = 0;
      return 0;
    }

}
void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{
  vector<string> nodes;
  vector<string> edges;
  
  // for(int i=0;i<nodes.size();i++)
  // 	cout << nodes[i]<< endl;
  // for(int i=0;i<edges.size();i++)
  // 	cout << edges[i]<< endl;
  int nb_perms[14][14];
  memset(nb_perms,0,sizeof(int[14][14]));
  int size_perms[14][14];
  memset(size_perms,0,sizeof(int[14][14]));

  char* files[14] = {"G0.ct","G1.ct","G2.ct","G3.ct","G4.ct","G5.ct","G6.ct","G7.ct","G8.ct","G9.ct","G10.ct","G11.ct","G12.ct","G13.ct",};
  for(int cur_file=0;cur_file<14;cur_file++){
    Graph3d* g = createGraph(files[cur_file],nodes,edges);
    treelet_spectrum** distrib = TreeletEnumerator::computeTreeletSpectrum(*g,nodes,edges);
    for(int i=0;i<cur_file;i++){
      vector<string> perms_nodes;
      int perm_node_size = 0;
      vector<string> perms_edges;
      int perm_edge_size = 0;
      if(!distrib[i]->empty()){
	for(treelet_spectrum::iterator it = distrib[i]->begin();it!=distrib[i]->end();it++){
	  perm_node_size = 0;
	  perm_edge_size = 0;
	  string code = it->first;
	  //Separation node/edge
	  stringstream ss_node;ss_node << "{";
	  stringstream ss_edge;ss_edge << "{";
	  bool is_node = true;
	  for(int c=0;c<code.length();c++)
	    if(code[c]!='_'){
	      if(is_node){
		ss_node << code[c];
		perm_node_size ++;
		if(c!=code.length()-1)
		  ss_node << ",";
	      }
	      else{
		ss_edge << code[c];
		perm_edge_size ++;
		if(c!=code.length()-3)
		  ss_edge << ",";
	      }
	      is_node =! is_node;
	    }
	  ss_node << "}";
	  perms_nodes.push_back(ss_node.str());
	  ss_edge << "}";
	  perms_edges.push_back(ss_edge.str());
	}
	cout << "const int TreeletEditDistance::G" << cur_file << "toG" << i ;
	cout << "[" << perms_nodes.size() << "][" << perm_node_size << "] = {";
	for(int cur_perm=0;cur_perm<perms_nodes.size();cur_perm++){
	  cout << perms_nodes[cur_perm];
	  if(cur_perm != perms_nodes.size()-1)
	    cout << ",";
	}
	cout << "}" << endl;
	cout << "TreeletEditDistance::edit_paths["<<cur_file<<"]["<<i<<"] = (void*)(TreeletEditDistance::G"<<cur_file<<"toG"<<i<<");" << endl;
	if(perm_edge_size > 0){
	  cout << "const int TreeletEditDistance::edgesG" << cur_file << "toG" << i ;
	  cout << "[" << perms_edges.size() << "][" << perm_edge_size << "] = {";
	  for(int cur_perm=0;cur_perm<perms_edges.size();cur_perm++){
	  cout << perms_edges[cur_perm];
	  if(cur_perm != perms_edges.size()-1)
	    cout << ",";
	  }
	  cout << "}" << endl;
	  cout << "TreeletEditDistance::edit_paths_edges["<<cur_file<<"]["<<i<<"] = (void*)(TreeletEditDistance::edgesG"<<cur_file<<"toG"<<i<<");" << endl;
	}
      }
      nb_perms[cur_file][i] = perms_nodes.size();
      size_perms[cur_file][i] = perm_node_size;
    }
  }
  cout << "nb_perms[14][14]={";
  for(int i =0;i<14;i++)
    {
      cout << "{";
      for(int j =0;j<14;j++){
	cout << nb_perms[i][j];
	if(j!=13)
	  cout << ",";
      }
      cout << "}," << endl;
    }
  cout << "}" << endl;

  cout << "size_perms[14][14]={";
  for(int i =0;i<14;i++)
    {
      cout << "{";
      for(int j =0;j<14;j++){
	cout << size_perms[i][j];
	if(j!=13)
	  cout << ",";
      }
      cout << "}," << endl;
    }
  cout << "}" << endl;

  return 0;
}

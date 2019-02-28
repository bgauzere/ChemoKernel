/*
 * @file MoleculeGraph.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *
 * @version 1.1.1 (2010-10-27)
 */

#include "MoleculeGraph.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <assert.h>
#include "separator.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

map<string, int> MoleculeGraph::table;
map<string, int> MoleculeGraph::shellTable;


MoleculeGraph::MoleculeGraph(const char* f, long id)
{
  if(loadGraph(f))
    _collection->SETVALUE("id", Long , (Long)id);
}

MoleculeGraph::~MoleculeGraph()
{
  if(_collection!=NULL)
    {

      delete[] _weights;
      delete[] _edges;

      delete[] _edgesIncidentNodes;

      for (int i=0; i<_collection->GETPARRAYSIZE("nodes", Collection); ++i)
	{
	  _nodes[i]->Delete();
	  delete _nodes[i];
	}
            
      delete[] _nodes;
      
      _collection->Delete();
      _graph->Delete();
      delete _graph;
      delete _collection;
    }
}

int MoleculeGraph::loadGraph(const char* f)
{
  this->filename = string(f);
  vector<char*> v;
  v = split(f, ".");
  if(!strcmp(v[v.size()-1],"ct"))
    {
      return loadGraphCT();
    }
  else
    {
      cerr << "Error : Can't recognize file format : " << this->filename << endl;
      _collection = NULL;
    }
  return 0;
}


int MoleculeGraph::loadGraphCT()
{
  
  // int nbBonds;
  int nbAtoms;
  ifstream file(this->filename.c_str(),ios::in);
  vector<char*> v;
  if (file.is_open()) 
    {		
      char * s = new char[255];
		
      file.getline(s, 255); // The first line is useless
      file.getline(s, 255); // Second line = NumberOfAtoms NumberOfBonds
      
      
      v = split(s, " ");
      
      nbAtoms = atoi(v[0]);
      nbBonds = atoi(v[1]);
      		
      _weights = new float[nbBonds];
      _graph = new Graph3d(false);
      _collection = new Collection;
      _nodes = new Collection*[nbAtoms];
      _edges = new char[nbBonds];
      _edgesIncidentNodes = new pair<int,int>[nbBonds];

      // Graph initialization
		
      _graph->New(nbAtoms, 0, 0, 0);
		
      // Creation of the nodes
		
      for(int i=0; i<nbAtoms; i++)
	{
	  file.getline(s,255); // s = x y z AtomLabel
	  v = split(s," ");
			
	  char index, c_atom;
	  string atom = v[3];
	  // if (atom.size() >= 2)
	  //   c_atom = atom.substr(0,2)[0];
	  // else
	  //   c_atom = atom.substr(0,1)[0];
	  c_atom = atom[0];
	  index = MoleculeGraph::table[atom];
	  
	  // cout << atom << endl;
	  // cout << (int)index << endl;
	  // cout << c_atom << endl;
	  // cout << "----" << endl;


	  _nodes[i] = new Collection;
	  _nodes[i]->SETVALUE("atom", Char, index);
	  _nodes[i]->SETVALUE("c_atom", Char, c_atom);
	  _graph->Add(i,i); // The i-th atom's index is at the i-th position of the nodes table.
	}
		
      // Creation of the edges
      for (int i=0; i<nbBonds; i++)
	{
	  file.getline(s,255); // s = Atom1 Atom2 BondType BondType
	  //vector<char*> v=split(s," ");
	  //Nicolo' Navarin: FIX reading
	  std::string s1=string(s);
	  int incident_node1, incident_node2;
	  istringstream buffer_node1(s1.substr(0,3));
	  buffer_node1 >> incident_node1;
	  istringstream buffer_node2(s1.substr(3,3));
	  buffer_node2 >> incident_node2;
	  
	  char label = (s1.substr(6,3).c_str())[0];
	  _edges[i] = label;
	  //cout<<"-------EDGE LABEL "<<v[2][0]<<endl;
	  _edgesIncidentNodes[i].first = min(incident_node1-1, incident_node2-1);
	  _edgesIncidentNodes[i].second = max(incident_node1-1, incident_node2-1);
	  
	 
	  _graph->Link(_edgesIncidentNodes[i].first, _edgesIncidentNodes[i].second,i,1.0f,false);
	  _weights[i] = 1.0f; // Needed by GraphToolBox
	}
		
      /* 
	 The "btrails-idx" index is used in order to avoid multiple bag of trails
	 computations of the same molecule. It is used in the file ShapeBagTrails.cpp
	 of the GraphToolBox library.
      */
		
      _collection->SETVALUE("btrails_idx", Long, -1);
		
      // We add the nodes, edges and weights tables and the graph into the collection.
      _collection->SETPARRAY("nodes", Collection, _nodes, nbAtoms);
      _collection->SETPOBJECT("graph", Graph3d, _graph);
      _collection->SETVALUE("sum_weight", Float, (float) nbBonds); // Needed by GraphToolBox
      _collection->SETARRAY("weights", Float, _weights, nbBonds);
      _collection->SETARRAY("edges", Char, _edges, nbBonds);
      
      file.close();
      return 1;
    }
  else
    {
      cerr << "Error : Can't open file : " << this->filename << endl;
      _collection = NULL;
      return 0;
    }
}


/** 
 * Les valeurs 123,124,125,126 sont utilisé par des séparateur.
 * Les séparateurs se trouvent dans separator.h
 */
void MoleculeGraph::initTable ()
{
  MoleculeGraph::table["H"] = 1;
  MoleculeGraph::table["He"] = 2;
  MoleculeGraph::table["Li"] = 3;
  MoleculeGraph::table["Be"] = 4;
  MoleculeGraph::table["B"] = 5;
  MoleculeGraph::table["C"] = 6;
  MoleculeGraph::table["N"] = 7;
  MoleculeGraph::table["O"] = 8;
  MoleculeGraph::table["F"] = 9;
  MoleculeGraph::table["Ne"] = 10;
  MoleculeGraph::table["Na"] = 11;
  MoleculeGraph::table["Mg"] = 12;
  MoleculeGraph::table["Al"] = 13;
  MoleculeGraph::table["Si"] = 14;
  MoleculeGraph::table["P"] = 15;
  MoleculeGraph::table["S"] = 16;
  MoleculeGraph::table["Cl"] = 17;
  MoleculeGraph::table["Ar"] = 18;
  MoleculeGraph::table["K"] = 19;
  MoleculeGraph::table["Ca"] = 20;
  MoleculeGraph::table["Sc"] = 21;
  MoleculeGraph::table["Ti"] = 22;
  MoleculeGraph::table["V"] = 23;
  MoleculeGraph::table["Cr"] = 24;
  MoleculeGraph::table["Mn"] = 25;
  MoleculeGraph::table["Fe"] = 26;
  MoleculeGraph::table["Co"] = 27;
  MoleculeGraph::table["Ni"] = 28;
  MoleculeGraph::table["Cu"] = 29;
  MoleculeGraph::table["Zn"] = 30;
  MoleculeGraph::table["Ga"] = 31;
  MoleculeGraph::table["Ge"] = 32;
  MoleculeGraph::table["As"] = 33;
  MoleculeGraph::table["Se"] = 34;
  MoleculeGraph::table["Br"] = 35;
  MoleculeGraph::table["Kr"] = 36;
  MoleculeGraph::table["Rb"] = 37;
  MoleculeGraph::table["Sr"] = 38;
  MoleculeGraph::table["Y"] = 39;
  MoleculeGraph::table["Zr"] = 40;
  MoleculeGraph::table["Nb"] = 41;
  MoleculeGraph::table["Mo"] = 42;
  MoleculeGraph::table["Tc"] = 43;
  MoleculeGraph::table["Ru"] = 44;
  MoleculeGraph::table["Rh"] = 45;
  MoleculeGraph::table["Pd"] = 46;
  MoleculeGraph::table["Ag"] = 47;
  MoleculeGraph::table["Cd"] = 48;
  MoleculeGraph::table["In"] = 49;
  MoleculeGraph::table["Sn"] = 50;
  MoleculeGraph::table["Sb"] = 51;
  MoleculeGraph::table["Te"] = 52;
  MoleculeGraph::table["I"] = 53;
  MoleculeGraph::table["Xe"] = 54;
  MoleculeGraph::table["Cs"] = 55;
  MoleculeGraph::table["Ba"] = 56;
  MoleculeGraph::table["La"] = 57;
  MoleculeGraph::table["Ce"] = 58;
  MoleculeGraph::table["Pr"] = 59;
  MoleculeGraph::table["Nd"] = 60;
  MoleculeGraph::table["Pm"] = 61;
  MoleculeGraph::table["Sm"] = 62;
  MoleculeGraph::table["Eu"] = 63;
  MoleculeGraph::table["Gd"] = 64;
  MoleculeGraph::table["Tb"] = 65;
  MoleculeGraph::table["Dy"] = 66;
  MoleculeGraph::table["Ho"] = 67;
  MoleculeGraph::table["Er"] = 68;
  MoleculeGraph::table["Tm"] = 69;
  MoleculeGraph::table["Yb"] = 70;
  MoleculeGraph::table["Lu"] = 71;
  MoleculeGraph::table["Hf"] = 72;
  MoleculeGraph::table["Ta"] = 73;
  MoleculeGraph::table["W"] = 74;
  MoleculeGraph::table["Re"] = 75;
  MoleculeGraph::table["Os"] = 76;
  MoleculeGraph::table["Ir"] = 77;
  MoleculeGraph::table["Pt"] = 78;
  MoleculeGraph::table["Au"] = 79;
  MoleculeGraph::table["Hg"] = 80;
  MoleculeGraph::table["Tl"] = 81;
  MoleculeGraph::table["Pb"] = 82;
  MoleculeGraph::table["Bi"] = 83;
  MoleculeGraph::table["Po"] = 84;
  MoleculeGraph::table["At"] = 85;
  MoleculeGraph::table["Rn"] = 86;
  MoleculeGraph::table["Fr"] = 87;
  MoleculeGraph::table["Ra"] = 88;
  MoleculeGraph::table["Ac"] = 89;
  MoleculeGraph::table["Th"] = 90;
  MoleculeGraph::table["Pa"] = 91;
  MoleculeGraph::table["U"] = 92;
  MoleculeGraph::table["Np"] = 93;
  MoleculeGraph::table["Pu"] = 94;
  MoleculeGraph::table["Am"] = 95;
  MoleculeGraph::table["Cm"] = 96;
  MoleculeGraph::table["Bk"] = 97;
  MoleculeGraph::table["Cf"] = 98;
  MoleculeGraph::table["Es"] = 99;
  MoleculeGraph::table["Fm"] = 100;
  MoleculeGraph::table["Md"] = 101;
  MoleculeGraph::table["No"] = 102;
  MoleculeGraph::table["Lr"] = 103;
  MoleculeGraph::table["Rf"] = 104;
  MoleculeGraph::table["Db"] = 105;
  MoleculeGraph::table["Sg"] = 106;
  MoleculeGraph::table["Bh"] = 107;
  MoleculeGraph::table["Hs"] = 108;
  MoleculeGraph::table["Mt"] = 109;
  MoleculeGraph::table["Ds"] = 110;
  MoleculeGraph::table["Rg"] = 111;
  MoleculeGraph::table["Cn"] = 112;
  MoleculeGraph::table["Uut"] = 113;
  MoleculeGraph::table["Uuq"] = 114;
  MoleculeGraph::table["Uup"] = 115;
  MoleculeGraph::table["Uuh"] = 116;
  MoleculeGraph::table["Uus"] = 117;
  MoleculeGraph::table["Uuo"] = 118;
  MoleculeGraph::table["D"] = 119; // Deuterium (isotope de H)


  MoleculeGraph::shellTable["H"] = 1;
  MoleculeGraph::shellTable["He"] = 1;
  MoleculeGraph::shellTable["Li"] = 2;
  MoleculeGraph::shellTable["Be"] = 2;
  MoleculeGraph::shellTable["B"] = 2;
  MoleculeGraph::shellTable["C"] = 2;
  MoleculeGraph::shellTable["N"] = 2;
  MoleculeGraph::shellTable["O"] = 2;
  MoleculeGraph::shellTable["F"] = 2;
  MoleculeGraph::shellTable["Ne"] = 2;
  MoleculeGraph::shellTable["Na"] = 3;
  MoleculeGraph::shellTable["Mg"] = 3;
  MoleculeGraph::shellTable["Al"] = 3;
  MoleculeGraph::shellTable["Si"] = 3;
  MoleculeGraph::shellTable["P"] = 3;
  MoleculeGraph::shellTable["S"] = 3;
  MoleculeGraph::shellTable["Cl"] = 3;
  MoleculeGraph::shellTable["Ar"] = 3;
  MoleculeGraph::shellTable["K"] = 4;
  MoleculeGraph::shellTable["Ca"] = 4;
  MoleculeGraph::shellTable["Sc"] = 4;
  MoleculeGraph::shellTable["Ti"] = 4;
  MoleculeGraph::shellTable["V"] = 4;
  MoleculeGraph::shellTable["Cr"] = 4;
  MoleculeGraph::shellTable["Mn"] = 4;
  MoleculeGraph::shellTable["Fe"] = 4;
  MoleculeGraph::shellTable["Co"] = 4;
  MoleculeGraph::shellTable["Ni"] = 4;
  MoleculeGraph::shellTable["Cu"] = 4;
  MoleculeGraph::shellTable["Zn"] = 4;
  MoleculeGraph::shellTable["Ga"] = 4;
  MoleculeGraph::shellTable["Ge"] = 4;
  MoleculeGraph::shellTable["As"] = 4;
  MoleculeGraph::shellTable["Se"] = 4;
  MoleculeGraph::shellTable["Br"] = 4;
  MoleculeGraph::shellTable["Kr"] = 4;
  MoleculeGraph::shellTable["Rb"] = 5;
  MoleculeGraph::shellTable["Sr"] = 5;
  MoleculeGraph::shellTable["Y"] = 5;
  MoleculeGraph::shellTable["Zr"] = 5;
  MoleculeGraph::shellTable["Nb"] = 5;
  MoleculeGraph::shellTable["Mo"] = 5;
  MoleculeGraph::shellTable["Tc"] = 5;
  MoleculeGraph::shellTable["Ru"] = 5;
  MoleculeGraph::shellTable["Rh"] = 5;
  MoleculeGraph::shellTable["Pd"] = 5;
  MoleculeGraph::shellTable["Ag"] = 5;
  MoleculeGraph::shellTable["Cd"] = 5;
  MoleculeGraph::shellTable["In"] = 5;
  MoleculeGraph::shellTable["Sn"] = 5;
  MoleculeGraph::shellTable["Sb"] = 5;
  MoleculeGraph::shellTable["Te"] = 5;
  MoleculeGraph::shellTable["I"] = 5;
  MoleculeGraph::shellTable["Xe"] = 5;
  MoleculeGraph::shellTable["Cs"] = 6;
  MoleculeGraph::shellTable["Ba"] = 6;
  MoleculeGraph::shellTable["La"] = 6;
  MoleculeGraph::shellTable["Ce"] = 6;
  MoleculeGraph::shellTable["Pr"] = 6;
  MoleculeGraph::shellTable["Nd"] = 6;
  MoleculeGraph::shellTable["Pm"] = 6;
  MoleculeGraph::shellTable["Sm"] = 6;
  MoleculeGraph::shellTable["Eu"] = 6;
  MoleculeGraph::shellTable["Gd"] = 6;
  MoleculeGraph::shellTable["Tb"] = 6;
  MoleculeGraph::shellTable["Dy"] = 6;
  MoleculeGraph::shellTable["Ho"] = 6;
  MoleculeGraph::shellTable["Er"] = 6;
  MoleculeGraph::shellTable["Tm"] = 6;
  MoleculeGraph::shellTable["Yb"] = 6;
  MoleculeGraph::shellTable["Lu"] = 6;
  MoleculeGraph::shellTable["Hf"] = 6;
  MoleculeGraph::shellTable["Ta"] = 6;
  MoleculeGraph::shellTable["W"] = 6;
  MoleculeGraph::shellTable["Re"] = 6;
  MoleculeGraph::shellTable["Os"] = 6;
  MoleculeGraph::shellTable["Ir"] = 6;
  MoleculeGraph::shellTable["Pt"] = 6;
  MoleculeGraph::shellTable["Au"] = 6;
  MoleculeGraph::shellTable["Hg"] = 6;
  MoleculeGraph::shellTable["Tl"] = 6;
  MoleculeGraph::shellTable["Pb"] = 6;
  MoleculeGraph::shellTable["Bi"] = 6;
  MoleculeGraph::shellTable["Po"] = 6;
  MoleculeGraph::shellTable["At"] = 6;
  MoleculeGraph::shellTable["Rn"] = 6;
  MoleculeGraph::shellTable["Fr"] = 7;
  MoleculeGraph::shellTable["Ra"] = 7;
  MoleculeGraph::shellTable["Ac"] = 7;
  MoleculeGraph::shellTable["Th"] = 7;
  MoleculeGraph::shellTable["Pa"] = 7;
  MoleculeGraph::shellTable["U"] = 7;
  MoleculeGraph::shellTable["Np"] = 7;
  MoleculeGraph::shellTable["Pu"] = 7;
  MoleculeGraph::shellTable["Am"] = 7;
  MoleculeGraph::shellTable["Cm"] = 7;
  MoleculeGraph::shellTable["Bk"] = 7;
  MoleculeGraph::shellTable["Cf"] = 7;
  MoleculeGraph::shellTable["Es"] = 7;
  MoleculeGraph::shellTable["Fm"] = 7;
  MoleculeGraph::shellTable["Md"] = 7;
  MoleculeGraph::shellTable["No"] = 7;
  MoleculeGraph::shellTable["Lr"] = 7;
  MoleculeGraph::shellTable["Rf"] = 7;
  MoleculeGraph::shellTable["Db"] = 7;
  MoleculeGraph::shellTable["Sg"] = 7;
  MoleculeGraph::shellTable["Bh"] = 7;
  MoleculeGraph::shellTable["Hs"] = 7;
  MoleculeGraph::shellTable["Mt"] = 7;
  MoleculeGraph::shellTable["Ds"] = 7;
  MoleculeGraph::shellTable["Rg"] = 7;
  MoleculeGraph::shellTable["Cn"] = 7;
  MoleculeGraph::shellTable["Uut"] = 7;
  MoleculeGraph::shellTable["Uuq"] = 7;
  MoleculeGraph::shellTable["Uup"] = 7;
  MoleculeGraph::shellTable["Uuh"] = 7;
  MoleculeGraph::shellTable["Uus"] = 7;
  MoleculeGraph::shellTable["Uuo"] = 7;
  MoleculeGraph::shellTable["D"] = 1; // Deuterium (isotope de H)
}

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

/*************************************************************/
/*Utils*/
/*************************************************************/

int MoleculeGraph::nbNeighbours(long i){
  int nb_neighbours = 0;
  for(GEdge * e =(*_graph)[i]->Neighbours();e!=NULL;e=e->Next(),nb_neighbours++)
    ;
  return nb_neighbours;
}

int MoleculeGraph::nbNeighbours(long i,Graph3d* g){
  int nb_neighbours = 0;
  for(GEdge * e =(*g)[i]->Neighbours();e!=NULL;e=e->Next(),nb_neighbours++)
    ;
  return nb_neighbours;
}

string  MoleculeGraph::indexToChemicalElement(char label){
  string chemical_element;
  for(map<string, int>::iterator it = MoleculeGraph::table.begin();
      it != MoleculeGraph::table.end(); it++)
    {
      if ((int) label == it->second)
	chemical_element = it->first;
    }
  return chemical_element;
}

string  MoleculeGraph::indexToChemicalElement(int label){
  string chemical_element;
  for(map<string, int>::iterator it = MoleculeGraph::table.begin();
      it != MoleculeGraph::table.end(); it++)
    {
      if ( label == it->second)
	chemical_element = it->first;
    }
  return chemical_element;
}

string MoleculeGraph::translateTreeletCode(string s){
  //By convention, we always start by a node
  stringstream cool_code;
  bool is_node =true;
  
  for(unsigned int i=0;i<s.size();i++)
    {
      if(s[i] != ' ' && s[i] != CHAR_SEP ){
	if(is_node)
	  cool_code << indexToChemicalElement(s[i]);
	else
	  cool_code << s[i];
	is_node = !is_node;
      }else
	if(s[i] == CHAR_SEP)
	  cool_code << LBL_SEP;
      // if(s[i] != ' '){

      // 	  cool_code.append(indexToChemicalElement(s[i]));
      // 	else
      // 	  cool_code.append(1,(char)s[i]);
	//cool_code.append("-");			 
      // }
      // else{
      // 	cool_code.append(" ");			 
      // 	is_node = true;
      // }
    }
  return cool_code.str();
  
}


string MoleculeGraph::translateTreeletCycleCode(string s){
  //By convention, we always start by a node
  stringstream cool_code;
  vector<char*> tokens = split(s.c_str(),LBL_SEP);
  for(unsigned int i=0;i<tokens.size();i++)
    {
      if(i!=0)
	cool_code << LBL_SEP;	
      bool is_node =true;
      for(unsigned int j=0;j<strlen(tokens[i]);j++)
	{
	  if(tokens[i][j] != ' '){
	    if(is_node)
	      cool_code << indexToChemicalElement(tokens[i][j]);
	    else
	      cool_code << tokens[i][j];
	    is_node = !is_node;
	  }
	}
    }
  return cool_code.str();
  
}

void MoleculeGraph::computeSpecialLabeledTreeletVector(treelet_spectrum ** specials,int nb_treelets){
  this->labeled_special = new CImg<double>(nb_treelets);
  assert(this->labeled_special->width() == nb_treelets);
  int cur_treelet = 0;
  for (int treelet_type=0;treelet_type<SIZE_SPECTRUM;treelet_type++){
    treelet_spectrum::iterator it = (*specials)[treelet_type].begin();
    for(;it != (*specials)[treelet_type].end();it ++){
      treelet_spectrum::iterator it_local = labeled_spectrum[treelet_type]->find(it->first);
      (*labeled_special)(cur_treelet) = (it_local == labeled_spectrum[treelet_type]->end())?0:it_local->second;
      cur_treelet ++;
    }
  }

  assert (cur_treelet == nb_treelets);
  assert (labeled_special != 0);
  long cast_labeled_special = (long)(this->labeled_special);
  _collection->SETVALUE("labeled_special", Long, (Long)(cast_labeled_special));
}
CImg<double> * MoleculeGraph::getSpecialLabeledTreeletSpectrum(){
  return labeled_special;
}

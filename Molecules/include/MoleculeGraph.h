/*
 * @file MoleculeGraph.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 * @author Benoit GAÜZÈRE <benoit.gauzere@ensicaen.fr>
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
*/

#ifndef __MOLECULE_GRAPH_H__
#define __MOLECULE_GRAPH_H__


#include <vector>
#include <map>
#include <set>
#include <string>
#include "pandore.h"
#include "CImg.h"
#include "TreeletEnumerator.h"
#include "TreeletEnumeratorAugmentedCycles.h"


// #define SIZE_SPECTRUM 6


/**
 * This class is a graph representation of a molecule.
 *
 * The informations about the molecule (atoms, bonds) are stored
 * as a graph into a Pandore collection (http://www.greyc.ensicaen.fr/~regis/Pandore/).
 */


class MoleculeGraph
{  	
  /*Attributes*/
  int nbBonds;
  pandore::Collection** _nodes; // The nodes table  
  std::pair<int,int> * _edgesIncidentNodes;
  char* _edges; // The bonds table
  float* _weights; // Used by GraphToolBox only
  pandore::Graph3d* _graph; // The graph representing the molecule
  pandore::Collection* _collection; // The Pandore collection
 
  /*Unlabeled Treelets (GbR 2011)*/
  double * _dgSpectrum; //The spectrum of Graphlets. _dgSpectrum[i] = nb graphlets_i
  int _sizeSpectrum; //The numbers of differents graphlets
  int nbPathsSizeK_recursive(pandore::GNode<pandore::Point3d> * n_start, int k, 
			     pandore::GNode<pandore::Point3d> * n_parent);
  int nbPathsSizeK(int k);
  
  /*Labeled Treelets*/
  treelet_spectrum ** labeled_spectrum;
  cimg_library::CImg<double> * labeled_special;
  // void enumerateKPaths(pandore::GNode<pandore::Point3d> * current_node, int size, 
  // 		       char * path_code, pandore::GNode<pandore::Point3d> * parent_node);
  // void enumerate3StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours);
  // void enumerate4StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours);
  // void enumerate5Star(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours);
  // void countGraphlet(char * code, int graphlet,double count_value = 1.0);
  
  /***************************************/
  /*Cycles */
  treelet_spectrum ** cycle_spectrum;
  //typedef std::vector<int> cycle; 
  class Substituent{
  public:
    int _position;
    int _hypergraph_node;
    std::string _code;
    Substituent():_position(-1),_hypergraph_node(-1),_code(""){
    }
  };
  class Cycle{
  public:
    bool * _edges_print; //Edges vector representation for cycles
    bool * _nodes_print; //Nodes vector representation for cycles
    std::string _code;
    int _size; //number of edges of cycle
    std::map<int, Substituent> _subs;//ensemble de subs, clé : hypergraph_node
    std::map< std::pair<int,int>, unsigned short> _angles;
    void computeAngles();
  };
  static bool cmp_cycle(Cycle c1, Cycle c2) { return (c1._size < c2._size);};
  std::vector<Cycle> _cycles;
  std::vector<Cycle> _relevantCycles; 
  std::vector<Cycle> _simpleCycles; 
  std::vector<std::string> _simpleCycleCodes; 

  pandore::Graph3d* _cyclicSystem;
  std::vector<std::string> _cyclicSystemBonds;

  int getAncestor(int node, int root, std::vector<int> *);
  void computeVismaraCycles();
  void recordCycle(std::vector <int> cycle_edges);
  std::string computeCycleCode(std::vector <int> cycle_edges);
  std::string computeCycleLinkCode(std::vector<int> cycle_edges);
  void computeRelevantCyclesGraphsBonds();
  void computeRelevantCycleGraph();
  std::vector<int> getConnexeEdgeSequence(int seed, std::vector<int> edges, bool * visited);
 
  /*Return if cyclic part. If true, cycles_list is filled with cycles containing node or edge*/
  bool isCyclicNode(int graph_node);
  bool isCyclicEdge(int graph_edge);
  /***************************************/
  /* Contracted Cycle Hypergraph*/
  treelet_spectrum ** cc_hypergraph_spectrum;
  pandore::Graph3d* _ccHypergraph;
  std::vector<std::string> _ccHypergraph_nodes_labels;
  std::vector<std::string> _ccHypergraph_edges_labels;
  std::vector<int> * graph_to_hypergraph; //graph nodes to hypergraph nodes mapping
  //Hyperedges storage:
  int *  edges_to_hyperedges;/*Correspondance between _graph edges et hyperedges :
			       hyperedges_index[i] is equals to -1 if
			       _edges[i] is not an hyperedge*/
  std::vector<int> hyperedges_to_edges;
  std::vector< std::pair<std::set<int>,std::set<int> > > hyperedges; // Vecteur de paires d'ensembles de sommets


  pandore::Graph3d* _contracted_ccHypergraph;
  std::vector<std::string> _contracted_ccHypergraph_nodes_labels;
  std::vector<std::string> _contracted_ccHypergraph_edges_labels;
  void computeContractedCycleHypergraph();
  void transformHypergraphToGraph();
  
  
  std::string computeContractedNodesCode(std::vector<int> set_of_nodes);
  std::string computeContractedNodesCode_rec(int current_node,bool * visited);
  int contractSetofNodes(std::vector<int> vec_u);
  treelet_spectrum ** cc_augmented_spectrum;
  std::vector<std::string> computeCycleSubs(unsigned int node_index);
  /*Misc*/
  int nbNeighbours(long i);	
  int nbNeighbours(long i,pandore::Graph3d* g);

  static std::map<std::string, int> table; // The correspondance table
  static std::map<std::string, int> shellTable; // The correspondance table of the shell of atoms


  std::vector< std::vector< std::pair <int, int> > > _twoconnexeComponents;/*lists of edges composing 
									     each 2-connexe components of the graph*/
  void computeTwoconnexeComponents();/*Compute _twoconnexeComponents */
  void biconnect(int v, int u); /*[Tarjan 1972]*/

  void computeTree(int atom,int numberOfTree);

  
public:
  std::string filename;

  MoleculeGraph() {};
  MoleculeGraph(const char* f, long = 1); //f: chemdraw Connection Table file
  ~MoleculeGraph();
  
  int loadGraph(const char* f);//use loadGraphCT  according to the extension of the file
  int loadGraphCT();// chemdraw Connection Table file, return load success

  pandore::Collection* getCollection() const
  {
    return _collection;
  }
	
  static void initTable (); //Initialization chemical element table
  //Compute user-friendly codes
  static std::string translateTreeletCode(std::string code);
  static std::string translateTreeletCycleCode(std::string code);
  static std::string indexToChemicalElement(char label);
  static std::string indexToChemicalElement(int label);
  
  int getCyclomaticNumber(){return nbBonds - nbAtoms()+1;};//Assume that we have 1 connected component
  int nbAtoms() const { return _graph->Size(); }
  int nbEdges() const { return nbBonds; }
  int getNbBiConnexe() { computeTwoconnexeComponents(); return  _twoconnexeComponents.size();}
  std::vector< std::vector< std::pair <int, int> > > getTwoconnexeComponents()
  {computeTwoconnexeComponents(); return _twoconnexeComponents;}
 
  /*********************************************************************/
  /*       Unlabeled Treelets, GBR 2011                                */
  /* Compute the Graphlet Spectrum for max size = 6                    */
  /* V1 adapted for alkanes (unlabeled, acyclic and max degree = 4)    */
  /*********************************************************************/
  void computeGraphletSpectrum();
  double * getSpectrum();/*The spectrum must be computed*/ 
  int getSpectrumSize();/*The spectrum must be computed*/

  /*********************************************************************/
  /*       Labeled Treelets, PRL 2011                                  */
  /* Compute the treelet Spectrum for max size = 6 and labeled treelets */
  /*********************************************************************/
  void computeLabeledTreeletSpectrum();
  treelet_spectrum ** getLabeledTreeletSpectrum();  /*The spectrum must be computed*/ 
  void normalizeLabeledSpectrum();  /*The spectrum must be computed*/         
                                    /*Normalize the treelet distribution by the numbers of occurences*/
  
  void computeSpecialLabeledTreeletVector(treelet_spectrum ** specials, int nb_treelets);
  cimg_library::CImg<double> * getSpecialLabeledTreeletSpectrum();  /*The special must be computed*/ 

  
  
  /*********************************************************************/
  /*       Cycles (ICPR 2012)                                          */
  /*                                                                   */
  /*********************************************************************/
  void computeRelevantCycles();
  // std::vector< Cycle >  getRelevantCycles();/*Returns lists of edges: each list is composed 
  // 					      of the indexes of the nodes of a relevant cycle*/
  
  void computeSimpleCycles(int k);//Relevant Cycles must be computed
  std::vector<std::string >  getSimpleCycleCodes();/*Returns lists of edges: each list is composed 
					    of the indexes of the nodes of a relevant cycle*/

  
  void computeCycleTreeletSpectrum();
  treelet_spectrum ** getCycleTreeletSpectrum();  /*The spectrum must be computed*/ 
  


  /*********************************************************************/
  /*       Contracted Cycles Hypergraph (GBR 2013)                     */
  /*                                                                   */
  /*********************************************************************/
  void computeContractedCycleGraphs();
  void computeContractedCycleSpectrum();
  treelet_spectrum ** getContractedCycleSpectrum();  /*The spectrum must be computed*/ 

  /*********************************************************************/
  /*       Augmented Cycles (To be released :))                        */
  /*                                                                   */
  /*********************************************************************/
  void computeAugmentedCycles(); // Compute substituent information
  void computeAugmentedCycleSpectrum();
  treelet_spectrum ** getAugmentedCycleSpectrum();  /*The spectrum must be computed*/ 
};

//TODO: A déplacer
std::vector<char*> split (const char* chaine, const char* sep);

#endif // __MOLECULE_GRAPH_H__

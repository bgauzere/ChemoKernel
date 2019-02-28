/**
 * @file TreeletEnumeratorAugmentedCycles.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __TREELETENUMERATORCYCLES_H__
#define __TREELETENUMERATORCYCLES_H__

#include <map>
#include <string>
// #include <gsl/gsl_sf_gamma.h>
#include "pandore.h"
//#include "MoleculeGraph.h"


#define SIZE_MAX 6

#define SIZE_SPECTRUM 14

typedef std::map<std::string, double,bool (*)(std::string, std::string)> treelet_spectrum;
typedef std::map< std::pair<int,int>, unsigned short> Cycle;
//TODO: Rajouter les edges dans TreeletToCollection
class TreeletEnumeratorAugmentedCycles
{
  /*Labeled Treelets*/
  /*Reverse labels order according to LBL_SEP split*/
  static std::string reverse_code(std::string code);
  static void enumerateKPaths(pandore::GNode<pandore::Point3d> * current_node, int size, 
			      std::string path_code, pandore::GNode<pandore::Point3d> * parent_node,
			      treelet_spectrum ** spectre, 
			      pandore::Graph3d& g, 
			      std::vector<std::string> nodes,
			      std::vector<std::string> edges,
			      bool * flags_visited,
			      std::vector<int> nodes_sequence,
			      std::vector<Cycle> cycles);
  
  static void enumerate3StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
  				    treelet_spectrum ** spectre, 
  				    pandore::Graph3d& g, 
  				    std::vector<std::string> nodes,
  				    std::vector<std::string> edges,
				    std::vector<int> nodes_sequence,
				    std::vector<Cycle> cycles);

  static void enumerate4StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
  				    treelet_spectrum ** spectre,
  				    pandore::Graph3d& g, 
  				    std::vector<std::string> nodes,
  				    std::vector<std::string> edges,
				    std::vector<int> nodes_sequence,
				    std::vector<Cycle> cycles);

  static void enumerate5Star(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
  			     treelet_spectrum ** spectre,
  			     pandore::Graph3d& g, 
  			     std::vector<std::string> nodes,
  			     std::vector<std::string> edges,
			     std::vector<int> nodes_sequence,
			     std::vector<Cycle> cycles);

  static void countGraphlet(std::string code, int graphlet,double count_value, treelet_spectrum ** spectre,
			    std::vector<int> nodes_sequence,std::vector<Cycle> cycles);

  static   std::string pairToChar(const std::pair< std::string,std::string > * pair_list, 
			 int nb_pair);
  
  static   int compPair(const void * l, const void * r);
  
  static std::pair<std::string, std::string> * pairSort(const std::pair<std::string, 
									std::string> * to_sort, 
							int nb_pair,
							std::vector<int> & nodes_sequence, 
							Cycle angles);
  struct substituent_to_compare{
    std::string _code;
    int _node;
    Cycle _angles;
    int _index;
  }; 
  static  bool compSubs(struct substituent_to_compare l , struct substituent_to_compare r);
 
public:
  static treelet_spectrum** computeAugmentedTreeletSpectrum(pandore::Graph3d& g, 
						   std::vector<std::string> nodes,
						   std::vector<std::string> edges,
						   std::vector<Cycle> cycles);
  // static int getNbAtoms(int treelet_type);

  // static pandore::Collection * TreeletToCollection(int treelet_type, std::string code);


};



#endif // __TREELETENUMERATORCYCLES_H__

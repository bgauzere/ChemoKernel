/**
 * @file TreeletEnumerator.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2012
 * 
 *  
 */

#ifndef __TREELETENUMERATOR_H__
#define __TREELETENUMERATOR_H__

#include <map>
#include <string>
//#include <gsl/gsl_sf_gamma.h>
#include "pandore.h"
// #include "MoleculeGraph.h"

#define SIZE_MAX 6

#define SIZE_SPECTRUM 14

typedef std::map<std::string, double,bool (*)(std::string, std::string)> treelet_spectrum;
//TODO: Rajouter les edges dans TreeletToCollection
class TreeletEnumerator
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
			      bool * flags_visited);
  
  static void enumerate3StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
  				    treelet_spectrum ** spectre, 
  				    pandore::Graph3d& g, 
  				    std::vector<std::string> nodes,
  				    std::vector<std::string> edges);
  static void enumerate4StarsGraphs(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
				    treelet_spectrum ** spectre,
  				    pandore::Graph3d& g, 
  				    std::vector<std::string> nodes,
  				    std::vector<std::string> edges);
  static void enumerate5Star(pandore::GNode<pandore::Point3d> * star_center, int nb_neighbours,
			     treelet_spectrum ** spectre,
			     pandore::Graph3d& g, 
			     std::vector<std::string> nodes,
			     std::vector<std::string> edges);
  static void countGraphlet(std::string code, int graphlet,double count_value, treelet_spectrum ** spectre);

  


public:
  static treelet_spectrum** computeTreeletSpectrum(pandore::Graph3d& g, 
						   std::vector<std::string> nodes,
						   std::vector<std::string> edges);
  static int getNbAtoms(int treelet_type);

  static pandore::Collection * TreeletToCollection(int treelet_type, std::string code);



};



#endif // __TREELETENUMERATOR_H__

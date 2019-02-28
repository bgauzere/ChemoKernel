/**
 * @file InterTreeletKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Fri Apr  6 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __INTERTREELETKERNEL_H__
#define __INTERTREELETKERNEL_H__

#include <pandore.h>
#include "CImg.h"
#include "GraphKernel.h"
#include "TreeletEnumerator.h"
#include "Kernel.h"
class InterTreeletKernel : public GraphKernel
{
  cimg_library::CImg<double> sim_treelets;
  Kernel* k_inter_graphs;
  int nb_treelets;
  std::vector< std::pair<int,int> > selected_pairs;
  void regul();
public:
  InterTreeletKernel(GraphKernel * k_inter_treelets, 
		     Kernel * k_inter_graphs, treelet_spectrum ** treelets, 
		     int nb_treelets );
  InterTreeletKernel(GraphKernel * k_inter_treelets, 
		     Kernel * k_inter_graphs, pandore::Collection ** treelets, 
		     int nb_treelets );
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  void setSimTreelets(cimg_library::CImg<double> sim);
  
  cimg_library::CImg<double> * getSimTreelets();

};

#endif // __INTERTREELETKERNEL_H__

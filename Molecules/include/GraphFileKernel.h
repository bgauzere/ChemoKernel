/*
 * @file GraphFileKernel.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @description : Computes the Gram Matrix from a file with the following format :
 * width height \n
 * value1 value2 ... valuen
 *
 * @version 1.0.0 
 */

#ifndef __GRAPH_FILE_KERNEL_H__
#define __GRAPH_FILE_KERNEL_H__

#include <pandore.h>
#include <deque>

#include "CImg.h"
#include "GraphKernel.h"

class GraphFileKernel : public GraphKernel
{
  cimg_library::CImg<double> Gram;
public:
  GraphFileKernel(char * gram_file);
  //check the id value in given Collection
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  
  void normalize();
  void RBF(double delta);
};

#endif // __GRAPH_FILE_KERNEL_H__

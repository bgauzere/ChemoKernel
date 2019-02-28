/**
 * @file CombinedKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb 29 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __COMBINEDKERNEL_H__
#define __COMBINEDKERNEL_H__

#include <pandore.h>
#include "GraphKernel.h"
#include <vector>

class CombinedKernel : public GraphKernel
{
  std::vector<GraphKernel *> _kernels;
  std::vector<double> _weights;

  /* 0 for an additive combination of kernel
   * 1 for an multiplicative combination of kernel
   */
  int _typeOfCombination; 

public:
  CombinedKernel(std::vector<GraphKernel *> kernels, std::vector<double> weights,int typeOfCombination=0);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);

};

#endif // __COMBINEDKERNEL_H__

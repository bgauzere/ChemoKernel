/**
 * @file GaussianKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GAUSSIANKERNEL_H__
#define __GAUSSIANKERNEL_H__

#include "Kernel.h"
#include "MoleculeGraph.h"

class GaussianKernel: public Kernel
{
  double _sigma;
public:
  GaussianKernel(double sigma = 1.0);
  // ~GaussianKernel();
  double operator()(double x, double y);
  
};

#endif // __GAUSSIANKERNEL_H__

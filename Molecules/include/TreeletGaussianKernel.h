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

#ifndef __TREELETGAUSSIANKERNEL_H__
#define __TREELETGAUSSIANKERNEL_H__

#include "MoleculeSpectrumKernel.h"
#include "MoleculeGraph.h"

class TreeletGaussianKernel: public MoleculeSpectrumKernel
{
  double _sigma;
public:
  TreeletGaussianKernel(double sigma = 1.0);
  // ~TreeletGaussianKernel();
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
  
};

#endif // __TREELETGAUSSIANKERNEL_H__

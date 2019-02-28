/**
 * @file InnerProductKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __INNERGAUSSIAN_H__
#define __INNERGAUSSIAN_H__
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
class InnerGaussianKernel : public MoleculeSpectrumKernel
{
  double _sigma;
public:
  InnerGaussianKernel(double sigma = 1.0);
  // ~GaussianKernel();
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
};

#endif // __INNERGAUSSIAN_H__

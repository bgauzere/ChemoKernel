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

#ifndef __POLYNOMIALKERNEL_H__
#define __POLYNOMIALKERNEL_H__
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
class PolynomialKernel : public MoleculeSpectrumKernel
{
  int _d;
public:
  PolynomialKernel(int d = 2);
  // ~GaussianKernel();
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
};

#endif // __POLYNOMIALKERNEL_H__

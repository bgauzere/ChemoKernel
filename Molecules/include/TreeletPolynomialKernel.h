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

#ifndef __TREELETPOLYNOMIALKERNEL_H__
#define __TREELETPOLYNOMIALKERNEL_H__
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
class TreeletPolynomialKernel : public MoleculeSpectrumKernel
{
  int _d;
public:
  TreeletPolynomialKernel(int d = 2);
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
};

#endif // __TREELETPOLYNOMIALKERNEL_H__

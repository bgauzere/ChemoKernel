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

#ifndef __TREELETINNERGAUSSIAN_H__
#define __TREELETINNERGAUSSIAN_H__
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
class TreeletInnerGaussianKernel : public MoleculeSpectrumKernel
{
  double _sigma;
public:
  TreeletInnerGaussianKernel(double sigma = 1.0);
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
};

#endif // __TREELETINNERGAUSSIAN_H__

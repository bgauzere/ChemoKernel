/**
 * @file RandomKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Mar 15 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Random Kernel for comparison
 * All necessary references.
 */


#ifndef __RANDOMKERNEL_H__
#define __RANDOMKERNEL_H__

#include "MoleculeSpectrumKernel.h"
#include "MoleculeGraph.h"

class RandomKernel : public MoleculeSpectrumKernel
{
public:
  double operator()(treelet_spectrum m1, treelet_spectrum m2,int treelet_type);

};

#endif // __RANDOMKERNEL_H__

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


#ifndef __TREELETRANDOMKERNEL_H__
#define __TREELETRANDOMKERNEL_H__

#include "MoleculeSpectrumKernel.h"
#include "MoleculeGraph.h"

class TreeletRandomKernel : public MoleculeSpectrumKernel
{
public:
  double operator()(treelet_spectrum m1, treelet_spectrum m2,int treelet_type);

};

#endif // __TREELETRANDOMKERNEL_H__

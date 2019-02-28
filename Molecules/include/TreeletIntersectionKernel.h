/**
 * @file IntersectionKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Mar 15 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Histogram Intersection Kernel as described by Barla et al. 2003, Histogram Intersection Kernel for Image Classification
 * All necessary references.
 */


#include "MoleculeSpectrumKernel.h"
#include "MoleculeGraph.h"

#ifndef __TREELETINTERSECTIONKERNEL_H__
#define __TREELETINTERSECTIONKERNEL_H__

class TreeletIntersectionKernel : public MoleculeSpectrumKernel
{
public:
  double operator()(treelet_spectrum m1, treelet_spectrum m2,int treelet_type);
};

#endif // __TREELETINTERSECTIONKERNEL_H__

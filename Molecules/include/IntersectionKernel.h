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


#include "Kernel.h"
#include "MoleculeGraph.h"

#ifndef __INTERSECTIONKERNEL_H__
#define __INTERSECTIONKERNEL_H__

class IntersectionKernel : public Kernel
{
public:
  double operator()(double x, double y);
};

#endif // __INTERSECTIONKERNEL_H__

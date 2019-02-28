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

#ifndef __INNERPRODUCTKERNEL_H__
#define __INNERPRODUCTKERNEL_H__
#include "MoleculeGraph.h"
#include "Kernel.h"

class InnerProductKernel : public Kernel
{
public:
  double operator()(double x, double y);
};

#endif // __INNERPRODUCTKERNEL_H__

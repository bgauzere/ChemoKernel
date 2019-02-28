/**
 * @file MinMaxKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * MinMax Kernel as described by Ralaivola 2005 
 * All necessary references.
 */

#ifndef __MINMAXKERNEL_H__
#define __MINMAXKERNEL_H__
#include <pandore.h>
#include "GraphKernel.h"

class MinMaxKernel: public GraphKernel
{
public:
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  
};

#endif // __MINMAXKERNEL_H__

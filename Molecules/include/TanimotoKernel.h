/**
 * @file TanimotoKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __TANIMOTOKERNEL_H__
#define __TANIMOTOKERNEL_H__
#include <pandore.h>
#include "GraphKernel.h"
#include  "MoleculeSpectrumKernel.h"

class TanimotoKernel : public GraphKernel
{
  MoleculeSpectrumKernel * k;
public:
  TanimotoKernel();
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  
};

#endif // __TANIMOTOKERNEL_H__

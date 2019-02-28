/**
 * @file TreeletJaccardKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Fri Feb  1 2013
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __TREELETJACCARDKERNEL_H__
#define __TREELETJACCARDKERNEL_H__

#include "MoleculeSpectrumKernel.h"
#include "MoleculeGraph.h"


class TreeletJaccardKernel: public MoleculeSpectrumKernel
{
public:
  TreeletJaccardKernel();
  double operator()(treelet_spectrum m1, treelet_spectrum m2, int treelet_type);
};

#endif // __TREELETJACCARDKERNEL_H__

/*
 * @file TreeletCycleKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Feb 20 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include "TreeletCycleKernel.h"

TreeletCycleKernel::TreeletCycleKernel(TreeletKernel * kt,CycleKernel * kc,
				       double weight_treelet,double weight_cycle){
  _kt = kt;
  _kc = kc;
  _weight_treelet = weight_treelet;
  _weight_cycle = weight_cycle;
}


double TreeletCycleKernel::operator() (pandore::Collection* c1, pandore::Collection* c2){
  return _weight_treelet*(*_kt)(c1,c2) + _weight_cycle*(*_kc)(c1,c2);
}

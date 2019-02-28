/**
 * @file TreeletCycleKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Feb 20 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __TREELETCYCLEKERNEL_H__
#define __TREELETCYCLEKERNEL_H__

#include <pandore.h>
#include <deque>
#include <vector>
#include <map>
#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculeGraph.h"
#include "TreeletKernel.h"
#include "CycleKernel.h"

//#define SIZE_SPECTRUM 14

class TreeletCycleKernel : public GraphKernel
{
  TreeletKernel * _kt;
  CycleKernel * _kc;
  bool _normalize;
  double _weight_treelet;
  double _weight_cycle;
public:
  TreeletCycleKernel(TreeletKernel * kt,CycleKernel * kc, 
		     double weight_treelet, double weight_cycle);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __TREELETCYCLEKERNEL_H__

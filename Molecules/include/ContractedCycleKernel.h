/**
 * @file ContractedCycleKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Dec 11 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __CONTRACTEDCYCLEKERNEL_H__
#define __CONTRACTEDCYCLEKERNEL_H__

#include <pandore.h>
#include <deque>
#include <vector>
#include <map>
#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
#include "TreeletType.h"

class ContractedCycleKernel : public GraphKernel
{
  MoleculeSpectrumKernel * k;
  double _sigma;  /*RBF tuning*/
  bool _normalize;
public:
 
  ContractedCycleKernel(KernelType kernel_type,double sigma = 1.0);
  void selectTreelets(std::vector<std::string> treelets[SIZE_SPECTRUM]);
  void weightTreelets(std::map<std::string, double,bool (*)(std::string, std::string)> treelets[SIZE_SPECTRUM]);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __CONTRACTEDCYCLEKERNEL_H__

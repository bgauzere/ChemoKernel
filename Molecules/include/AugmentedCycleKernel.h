/*
 * @file AugmentedCycleKernel.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @TODO: Destructeur
 * @version 1.0.0 (2010-10-27)
 */

#ifndef __AUGMENTED_CYCLE_KERNEL_H__
#define __AUGMENTED_CYCLE_KERNEL_H__

#include <pandore.h>
#include <deque>
#include <vector>
#include <map>
#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculeGraph.h"
#include "MoleculeSpectrumKernel.h"
#include "TreeletType.h"

class AugmentedCycleKernel : public GraphKernel
{
  MoleculeSpectrumKernel * k;
  double _sigma;  /*RBF tuning*/
  bool _normalize;
public:

  AugmentedCycleKernel(KernelType kernel_type,double sigma = 1.0);
  void selectTreelets(std::vector<std::string> treelets[SIZE_SPECTRUM]);
  void weightTreelets(std::map<std::string, double,bool (*)(std::string, std::string)> treelets[SIZE_SPECTRUM]);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __TREELET_KERNEL_H__

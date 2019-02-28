/**
 * @file ContractedCycleKernel.cpp
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.0.0 (2010-10-27)
 */

#include <algorithm>
#include <cmath>
#include <cassert>
#include "ContractedCycleKernel.h"
#include "MoleculeGraph.h"
#include "TreeletGaussianKernel.h"
#include "TreeletInnerProductKernel.h"
#include "string_utils.h"
#include "TreeletBinaryKernel.h"
#include "TreeletIntersectionKernel.h"
#include "TreeletRandomKernel.h"
#include "TreeletCompleteGaussianKernel.h"
#include "TreeletInnerGaussianKernel.h"
#include "TreeletPolynomialKernel.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

ContractedCycleKernel::ContractedCycleKernel(KernelType kernel_type, double sigma){
  switch(kernel_type)
    {
    case IntersectionKernelType:
      k = new TreeletIntersectionKernel();
      break;
    case GaussianKernelType:
      k = new TreeletGaussianKernel(sigma);
      break;
    case InnerProductKernelType:
      k = new TreeletInnerProductKernel();
      break;
    case BinaryKernelType:
      k = new TreeletBinaryKernel();
      break;
    case RandomKernelType:
      k = new TreeletRandomKernel();
      break;
    case CompleteGaussianKernelType:
      k = new TreeletCompleteGaussianKernel(sigma);
      break;
    case InnerGaussianKernelType:
      k = new TreeletInnerGaussianKernel(sigma);
      break;
    case PolynomialKernelType:
      k = new TreeletPolynomialKernel(sigma);
      break;
    default:
      //devrait pas arriver
      break;
    }
}
void ContractedCycleKernel::selectTreelets(vector<string> treelets[SIZE_SPECTRUM])
{
  k->selectTreelets(treelets);
}
void ContractedCycleKernel::weightTreelets(map<string, double,bool (*)(string, string)> treelets[SIZE_SPECTRUM])
{
  k->weightTreelets(treelets);
}

double ContractedCycleKernel::operator() (Collection* c1, Collection* c2)
{
  double scalar = 0.0;
  treelet_spectrum ** spec_1 = (treelet_spectrum **) c1->GETVALUE("cc_hypergraph_spectrum",Long);
  treelet_spectrum ** spec_2 = (treelet_spectrum **) c2->GETVALUE("cc_hypergraph_spectrum",Long);
  
  for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
    {
      scalar += (*k)(spec_1[i][0],spec_2[i][0],i);
    }

  return scalar;  
}

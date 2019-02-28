/**
 * @file TreeletKernel.cpp
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.0.0 (2010-10-27)
 */

#include <algorithm>
#include <cmath>
#include <cassert>
#include "TreeletKernel.h"
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
#include "TreeletJaccardKernel.h"


#include <iostream>

using namespace std;
using namespace cimg_library;
using namespace pandore;
 
TreeletKernel::TreeletKernel(KernelType kernel_type, double sigma){
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
    case JaccardKernelType:
      k = new TreeletJaccardKernel();
      break;

    default:
      //devrait pas arriver
      break;
    }

  //k = new GaussianKernel(sigma);
  //k = new InnerProductKernel();
  //k = new BinaryKernel();

  //k = new RandomKernel();
  //_sigma = sigma;
}

void TreeletKernel::selectTreelets(vector<string> treelets[SIZE_SPECTRUM])
{
  k->selectTreelets(treelets);
}
void TreeletKernel::weightTreelets(map<string, double,bool (*)(string, string)> treelets[SIZE_SPECTRUM])
{
  k->weightTreelets(treelets);
}

double TreeletKernel::operator() (Collection* c1, Collection* c2)
{
  double scalar = 0.0;
  treelet_spectrum ** spec_1 = (treelet_spectrum **) c1->GETVALUE("labeled_spectrum",Long);
  treelet_spectrum ** spec_2 = (treelet_spectrum **) c2->GETVALUE("labeled_spectrum",Long);
  
  for(int i = 0; i < SIZE_SPECTRUM // 1
	; i++) //Parcours des 13 Treelets
    {
      scalar += (*k)(spec_1[i][0],spec_2[i][0],i);
      assert(!std::isnan(scalar));
      // treelet_spectrum::iterator it1 = spec_1[i]->begin();
      // while(it1 != spec_1[i]->end())
      // 	{
      // 	  if(spec_2[i]->find(it1->first) != spec_2[i]->end())
      // 	    {
      // 	      double diff = fabs(spec_2[i]->find(it1->first)->second - it1->second);
      // 	      scalar += exp(-(diff * diff)/(2*_sigma*_sigma));
      // 	    }
      // 	  it1++;
      // 	}      
    }

  /*Nombre de treelets*/
  // if(_normalize)
  //   {
  //     double N1 = 0;
  //     double N2 = 0;
  //     for(int i = 0; i < 13; i++)
  // 	{
  // 	  for(treelet_spectrum::iterator it1 = spec_1[i]->begin();it1!= spec_1[i]->end();it1++)
  // 	    { 
  // 	      N1 += it1->second;
  // 	    }
  // 	  for(treelet_spectrum::iterator it2 = spec_2[i]->begin();it2!= spec_2[i]->end();it2++)
  // 	    { 
  // 	      N2 += it2->second;
  // 	    }
  // 	}
  //     scalar /= sqrt(N1*N2);
  //   }
  
  //Nombre d'atomes
  // if(_normalize)
  //   {
  //     unsigned int n1 = c1->GETPARRAYSIZE("nodes", Collection);
  //     unsigned int n2 = c2->GETPARRAYSIZE("nodes", Collection);
  //     scalar /= sqrt(n1*n2);
  //   }

  return scalar;  
}

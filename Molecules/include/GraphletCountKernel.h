/*
 * @file GraphletCountKernel.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * 
 * @version 1.0.0 (2010-10-27)
 */

#ifndef __GRAPHLET_COUNT_KERNEL_H__
#define __GRAPHLET_COUNT_KERNEL_H__

#include <pandore.h>
#include <deque>

#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculeGraph.h"

//#define SIZE_SPECTRUM 14

class GraphletCountKernel : public GraphKernel
{
  double _sigma;  /*RBF tuning*/
  double _mu;     /* Power weight on graphlet weighting*/
  double * _graphlet_coeff /*Graphlet weighting*/;
public:
  GraphletCountKernel(double sigma = 1.0,double mu = 1.0);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
  void setCoeffs(double * coeffs, int size);
};

#endif // __GRAPHLET_COUNT_KERNEL_H__

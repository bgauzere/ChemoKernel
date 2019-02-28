/**
 * @file KernelRidge.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "KernelRidge.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
using namespace cimg_library;
using namespace pandore;

double KernelRidge::operator() (Collection* col)
{
  unsigned int n = _dataset->size();
  CImg<double> K (1,n);
  CImg<double> y (n,1);
	
  for (unsigned int i=0; i<n; ++i)
    {
      double k = (*_kg)(_dataset->getCollection(i), col);
      
      double a = (*_kg) (col, col);
      double b = (*_kg) (_dataset->getCollection(i),_dataset->getCollection(i));
      k /= sqrt(a*b);	
     
      K(0,i) = k;
      y(i,0) = _dataset->getParameter(i);
    }
	
  CImg<double> tmp = _dataset->getGramMatrix(true) + CImg<double>(n,n).identity_matrix()*_lambda;
  CImg<double> tmp2 = tmp.invert();
  CImg<double> res = y*tmp2*K;
	
  return res(0,0);
}

double KernelRidge::optimalParameter (Collection* col, double p)
{
  double lambda = 1.0;
  double max_lambda = 400.0;
  double epsilon = 5.0;
	
  _lambda = lambda;
  double val = (*this)(col);
  double delta = abs(val-p);
  double optimal_lambda = lambda;
	
  lambda += epsilon;
	
  while (lambda <= max_lambda)
    {
      _lambda = lambda;
      val = (*this)(col);
      double d = abs(val-p);
		
      if (d < delta)
	{
	  delta = d;
	  optimal_lambda = lambda;
	}
      lambda += epsilon;
    }
	
  return optimal_lambda;
}

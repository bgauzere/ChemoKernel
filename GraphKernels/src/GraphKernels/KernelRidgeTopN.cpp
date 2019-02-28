/**
 * @file KernelRidgeTopN.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "KernelRidgeTopN.h"
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
using namespace std;
using namespace cimg_library;
using namespace pandore;

bool compare ( pair<double, int > * a, pair<double, int > * b)
{
  return ( a->first > b->first);
}


KernelRidgeTopN::KernelRidgeTopN (GraphKernel* kg, Dataset* dataset, double lambda, int nb_neighbours)
  {
    _kg = kg;
    _dataset = dataset;
    _lambda = lambda;
    _nb_neighbours = nb_neighbours;
  }
  
double KernelRidgeTopN::operator() (Collection* col)
{
  unsigned int n = _dataset->size();
	
  CImg<double> K (1,_nb_neighbours);
  CImg<double> y (_nb_neighbours,1);
  
  vector< pair<double, int > * > neighbours_id;
  for (int i=0; i < n; i++)
    {
      neighbours_id.push_back(new pair<double, int>((*_kg)(_dataset->getCollection(i), col), i));
    }
  
  
  sort (neighbours_id.begin(), neighbours_id.end(), compare);
  // for (unsigned int i=0; i<neighbours_id.size(); ++i)
  //   {
  //     cout << neighbours_id[i]->second << " : " << neighbours_id[i]->first << endl;
  //   }
  
  for (unsigned int i=0; i<_nb_neighbours; ++i)
    {
      double k = neighbours_id[i]->first; //(*_kg)(_dataset->getCollection(neighbours_id[i]), col);
      
      double a = (*_kg) (col, col);
      double b = (*_kg) (_dataset->getCollection(neighbours_id[i]->second),
			 _dataset->getCollection(neighbours_id[i]->second));
      k /= sqrt(a*b);
      
      K(0,i) = k;
      y(i,0) = _dataset->getParameter(neighbours_id[i]->second);
    }
  
  CImg<double> gram = _dataset->getGramMatrix(true);
  CImg<double> gram_top (_nb_neighbours,_nb_neighbours);
  for(int i=0;i < _nb_neighbours;i++)
    {
      for(int j=0;j < _nb_neighbours;j++)
	{
	  gram_top(i,j) = gram(neighbours_id[i]->second,neighbours_id[j]->second);
	}
    }
  
  CImg<double> tmp = gram_top + CImg<double>(_nb_neighbours,_nb_neighbours).identity_matrix()*_lambda;
  // //Verification de tmp
  // CImg<double> eigvals;
  // CImg<double> eigvects;
  // tmp.symmetric_eigen(eigvals,eigvects);
  // bool isPD = true;
  // for (unsigned int i=0; i<_nb_neighbours; ++i)
  //   {
  //     cout << eigvals[i];
  //     // cout << eigvects[i] << endl;
  //     //   if(eigvals[i] < 0)
  //     // 	return false;
  //     // 
  //     if(eigvals[i] < 0)
  // 	cout << "Aïe aïe aïe";
  //     cout << endl;
  //   }

  CImg<double> tmp2 = tmp.invert();
  CImg<double> res = y*tmp2*K;
	
  return res(0,0);
}

void KernelRidgeTopN::setNbNeighbours(int nb_neighbours){
  _nb_neighbours = nb_neighbours;
}
double KernelRidgeTopN::optimalParameter (Collection* col, double p)
{
  double lambda = 0.001;
  double max_lambda = 20;
  double epsilon = 0.001;
	
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

/**
 * @file LaplacianKernelOriginal.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include <algorithm>
#include "LaplacianKernelOriginal.h"
#include "GraphEditDistance.h"
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <iomanip>
#include <sys/times.h>
#include <sys/time.h>

#include <unistd.h>

using namespace std;
using namespace cimg_library;
using namespace pandore;


LaplacianKernelOriginal::LaplacianKernelOriginal(const LaplacianKernelOriginal * clone){
  edit = clone->edit;
  sigma = clone->sigma;
  regularization = clone->regularization;
  lambda = clone->lambda;
  K = clone->K;
  W = clone->W;  
  D = clone->D;  
  L = clone->L;
  U_K = clone->U_K;
  _collections = clone->_collections;
}

LaplacianKernelOriginal::LaplacianKernelOriginal (GraphEditDistance* edit, const Dataset & trainset, const Dataset & testset, double sigma, int regularization, double lambda) : edit(edit), sigma(sigma), regularization(regularization), lambda(lambda)
{	
  unsigned int N = trainset.size();
  unsigned int M = testset.size();
	
  for (unsigned int i=0; i<N; ++i)
    _collections.push_back(trainset.getCollection(i));
	
  for (unsigned int i=0; i<M; ++i)
    _collections.push_back(testset.getCollection(i));
	
  computeKernel();
}

LaplacianKernelOriginal::LaplacianKernelOriginal (GraphEditDistance* edit, const Dataset & trainset, double sigma, int regularization,  double lambda) : edit(edit), sigma(sigma), regularization(regularization), lambda(lambda)
{	
  unsigned int N = trainset.size();
	
  for (unsigned int i=0; i<N; ++i)
    _collections.push_back(trainset.getCollection(i));
	
  computeKernel();
}

double LaplacianKernelOriginal::computeWeight(pandore::Collection* m1, pandore::Collection* m2){
  double dist = (*edit) (m1, m2);
  return exp(-dist/(sigma*sigma));
}

void LaplacianKernelOriginal::computeKernel ()
{	
  //Calcul complet
  unsigned int N = _collections.size();
	
  W.resize(N,N); // The weights matrix
  D.resize(N,N); // The degree matrix
  L.resize (N,N); // The graph Laplacian matrix
  U_K.resize(N, N);

  K.resize(N, N);
	
  // Computing the weights and degree matrix

  for (unsigned int i = 0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
	{
	  D(i,j) = 0;
	  if (i != j)
	    {
	      W(i,j) = computeWeight(_collections[i], _collections[j]);
	    }
	  else
	    W(i,j)=0;
	}
  
  for (unsigned int i=0; i<N; ++i)
    {
      double d = 0;
      for (unsigned int j=0; j<N; ++j)
	d += W(i,j);
		
      D(i,i) = d;
    }
  // Computing the graph Laplacian matrix
  //Normalized Version
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; j++)
      {
  	if (i==j)
  	  L(i,j)=1;
  	else if ((i!=j) && (W(i,j) != 0))
  	  L(i,j)=-W(i,j)/sqrt(D(i,i)*D(j,j));
  	else 
  	  L(i,j)=0;
      }

  // Computing the Laplacian Moore-Penrose pseudoinverse
  CImg<double> eigvals;
  CImg<double> eigvects;

  L.symmetric_eigen(eigvals,eigvects);

  K=(double) 0;
	
  for (unsigned int i=0; i<N; ++i)
    {
      //cout << eigvals[i] << endl;
      double r;
      switch (regularization)
  	{
  	case 0:
  	  r = eigvals[i]; // No regularization
  	  break;
  	case 1:
  	  r = 1 +lambda*eigvals[i]; // Regularized Laplacian
  	  break;
  	case 2:
  	  r = exp(eigvals[i]*lambda*lambda/2); // Diffusion Process
  	  break;
  	case 3:
  	  r = 1/(cos(eigvals[i]*M_PI/4)); // Inverse cosine
  	  break;
  	case 4:
  	  r = 1/(lambda-eigvals[i]); // One-step Random Walk
  	  break;
  	}
      assert(fabs(r) > 0.00001);
      CImg<double> tmp (N,N);
      for (unsigned int j=0; j<N; ++j)
  	tmp(j,j)= eigvects(i,j)*eigvects(i,j);
      for (unsigned int j=0; j<N; ++j)
  	for (unsigned int k=j+1; k<N; ++k)
  	  {
  	    tmp(j,k) = eigvects(i,j)*eigvects(i,k);
  	    tmp(k,j)=tmp(j,k);
  	  }
      K += (1/r)*tmp;
      //U_K += (1/r)*tmp;
    }

  
  //Normalized Laplacian
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      K(i,j) = sqrt(D(i,i))*K(i,j)*sqrt(D(j,j));
  
  normalizeK();
  
}

void LaplacianKernelOriginal::normalizeK(){
  //  Normalization
  unsigned int N  = _collections.size();
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      if (i != j)
  	K(i,j) /= sqrt(K(i,i)*K(j,j));  

  for (unsigned int i=0; i<N; ++i)
    K(i,i)=1;
}

 double LaplacianKernelOriginal::operator() (Collection* c1, Collection* c2)
{	
  int i=-1;
  int j=-1;
 	
  for (unsigned int k=0; k<_collections.size(); ++k)
    if (_collections[k] == c1)
      {
	i=k;
	break;
      }
  for (unsigned int k=0; k<_collections.size(); ++k)
    if (_collections[k] == c2)
      {
	j=k;
	break;
      }
	
  bool compute = false;

  if (i == -1)
    {
      i=_collections.size();
      _collections.push_back(c1);
      compute = true;
    }
	
  if (j == -1)
    {
      j=_collections.size();
      _collections.push_back(c2);
      compute = true;
    }
	
  if (compute) //The kernel needs to be recomputed
    {  
      clock_t start,end;
      struct timeval s_start, s_end;
      start = gettimeofday(&s_start,NULL);

      end = gettimeofday(&s_end,NULL);//times(&s_end);
      computeKernel();
      long diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Slow Add (micros) : " << diff << endl;
    }
  return K(i,j);
}

void LaplacianKernelOriginal::printWeightMatrix(){
  unsigned int N = _collections.size();
  cout.precision(DBL_DIG);
  for (unsigned int i = 0; i<N; ++i)
    {
      for (unsigned int j=0; j<N; ++j)
	{
	  double dist = (*edit) (_collections[i], _collections[j]);
	  cout << exp(-dist/(sigma*sigma)) << " ";
	  
	}
      cout << endl;
    }
}

void LaplacianKernelOriginal::printMatrix(CImg<double> m){
  if(PRINT_MATRIX)
    {
      for (int i=0; i< m.width(); ++i)
	{
	  for (int j=0; j< m.height(); ++j)
	    {
	      cout.precision(DBL_DIG);
	      cout <<  setw(20) << setiosflags(ios::left) << m(i,j)  << "\t";
	    }
	  cout << endl;
	}
    }
}

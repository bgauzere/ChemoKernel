/**
 * @file KernelPCA.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "KernelPCA.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

KernelPCA::KernelPCA (GraphKernel* kgraph, const Dataset & dataset) : kgraph(kgraph)
{	
  for (unsigned int i=0; i<dataset.size(); ++i)
    cols.push_back(dataset.getCollection(i));
	
  computeGramMatrix();
}

void KernelPCA::computeGramMatrix ()
{
  unsigned int N = cols.size();
	
  // We compute the Gram matrix
	
  k.resize(N,N);
	
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      k(i,j) = (*kgraph) (cols[i], cols[j]);
	
  // We compute the centered Gram matrix
	
  kCentered.resize(N,N);
	
  double tmp1 = 0.0;
	
  for (unsigned int r=0; r<N; ++r)
    for (unsigned int s=0; s<N; ++s)
      tmp1 += k(r,s);
	
  tmp1/=N*N;
	
  for (unsigned int i=0; i<N; ++i)
    {
      for (unsigned int j=0; j<N; ++j)
	{
	  kCentered(i,j) = k(i,j)+tmp1;
			
	  double tmp = 0.0;
	  for (unsigned int r=0; r<N; ++r)
	    tmp += k(i,r)+k(r,j);
			
	  tmp/=N;
			
	  kCentered(i,j) -= tmp;
	}
    }
	
  // We compute the eigenvalues and eigenvectors of the centered gram matrix
	
  kCentered.symmetric_eigen(eigvals,eigvects);
	
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
	if (eigvals(i) < 0.001)
	  eigvects(i,j) = 0;
	else
	  eigvects(i,j) = eigvects(i,j)/sqrt(eigvals(i));
      }
}

void KernelPCA::addToTrainingSet (Collection* col, bool compute)
{	
  cols.push_back(col);
	
  if (compute)
    {
      unsigned int N = cols.size();
		
      // We resize the Gram matrix
		
      k.resize(N,N);
		
      for (unsigned int i=0; i<N; ++i)
	for (unsigned int j=0; j<N-1; ++j)
	  k(i,j)=k(i,j+1);
		
      for (unsigned int i=0; i<N-1; ++i)
	for (unsigned int j=0; j<N; ++j)
	  k(i,j)=k(i+1,j);
		
      for (unsigned int i=0; i<N; ++i)
	{
	  k(i,N-1)=(*kgraph)(cols[i],cols[N-1]);
	  k(N-1,i)=k(i,N-1);
	}
		
      // We compute the centered Gram matrix
		
      kCentered.resize(N,N);
		
      double tmp1 = 0.0;
		
      for (unsigned int r=0; r<N; ++r)
	for (unsigned int s=0; s<N; ++s)
	  tmp1 += k(r,s);
		
      tmp1/=N*N;
		
      for (unsigned int i=0; i<N; ++i)
	{
	  for (unsigned int j=0; j<N; ++j)
	    {
	      kCentered(i,j) = k(i,j)+tmp1;
				
	      double tmp = 0.0;
	      for (unsigned int r=0; r<N; ++r)
		tmp += k(i,r)+k(r,j);
				
	      tmp/=N;
				
	      kCentered(i,j) -= tmp;
	    }
	}
		
      // We compute the eigenvalues and eigenvectors
		
      kCentered.symmetric_eigen(eigvals,eigvects);	
    }
}

double** KernelPCA::projectTrainingSet (unsigned int p) const
{
  unsigned int N = cols.size();
	
  double** pc = new double*[N];
	
  for (unsigned int i=0; i<N; ++i)
    {
      pc[i] = new double[p];
		
      for (unsigned int q=0; q<p; ++q)
	pc[i][q]=0.0;
		
      for (unsigned int j=0; j<N; ++j)
	{
	  double x = kCentered(i,j);
		  
	  for (unsigned int q=0; q<p; ++q)
	    pc[i][q] += eigvects(q,j)*x;
	}
    }	
	
  return pc;
}

double* KernelPCA::project (Collection* col, unsigned int p)
{
  double* pc = new double[p];
	
  unsigned int N = cols.size();
	
  double* buffer = new double[N];
	
  for (unsigned int i=0; i<N; ++i)
    buffer[i] = (*kgraph) (cols[i], col);
	
  unsigned int m = p;
	
  if (cols.size() < p)
    m = cols.size();
	
  for (unsigned int q=0; q<m; ++q)
    pc[q] = 0.0;
	
  double tmp1 = 0.0;
  for (unsigned int r=0; r<N; ++r)
    for (unsigned int s=0; s<N; ++s)
      tmp1 += k(r,s);
  tmp1/=N*N;
		
  for (unsigned int i=0; i<N; ++i)
    {
      double x = buffer[i]+tmp1;
			
      double tmp = 0.0;
      for (unsigned int r=0; r<N; ++r)
	tmp += k(i,r)+buffer[r];
      tmp/=N;
			
      x -= tmp;
			
      for (unsigned int q=0; q<m; ++q)
	pc[q] += eigvects(q,i)*x;
    }
	
  delete[] buffer;
	
  return pc;
}

double KernelPCA::reconstructionError (Collection* col, unsigned int p)
{
  unsigned int N = cols.size();
	
  // We project the graph onto the p principal eigenvectors
	
  double* pc = project(col, p);
	
  double ps = (*kgraph)(col, col);
	
  double tmp = 0.0;
	
  for (unsigned int i=0; i<N; ++i)
    tmp += (*kgraph)(col, cols[i]);
	
  ps -= (2/N)*tmp;
	
  tmp = 0.0;
	
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      tmp += k(i,j);
	
  ps += (1/(N*N))*tmp;
	
  tmp = 0.0;
	
  unsigned int m = p;
	
  if (cols.size() < p)
    m = cols.size();
	
  for (unsigned int q=0; q<m; ++q)
    {
      double x = pc[q];
		
      tmp += x*x;
    }
	
  delete[] pc;
	
  return (ps - tmp);
}

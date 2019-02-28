/**
 * @file Dataset.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include <cfloat>
#include <cassert>
#include "Dataset.h"
#include <fstream>
#include <math.h>

#include <iostream>
#include <sstream>

using namespace std;
using namespace cimg_library;
using namespace pandore;

Dataset::Dataset(const Dataset& d)
{
  _collections = d._collections;
  _parameters = d._parameters;
  _gram = d._gram;
}

void Dataset::computeGramMatrix(GraphKernel* kg, bool normalize)
{
  unsigned int n = _collections.size();
	
  _gram.resize (n,n,1,1);
	
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<=i; ++j)
      {
	//cout<<"(i,j)=("<<i<<","<<j<<")"<<endl;
	double x = (*kg) (_collections[i], _collections[j]);
			
	_gram(i,j) = x;
	_gram(j,i) = x;
      }
	
  if (normalize) // We normalize the Gram matrix
    {
      for (unsigned int i=0; i<n; ++i)
	for (unsigned int j=0; j<n; ++j)
	  if (i != j)
	    _gram(i,j) /= sqrt(_gram(i,i)*_gram(j,j));
		
      for (unsigned int i=0; i<n; ++i)
	_gram (i,i)=1;
    }
}

void Dataset::showIdenticalLine()
{
  unsigned int n = _collections.size();
  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int i2=i+1; i2<n; ++i2)
	{
	  unsigned int j=0;
	  while(j<n && _gram(i,j)==_gram(i2,j))
	    j++;
	  if(j==n)
	    cout<<"Line "<<i<<" and "<<i2<<" are identical"<<endl;
	}
    } 
}

CImg<double> Dataset::getGramMatrix (bool normalize) const
{
  if (!normalize)
    return _gram;
	
  // We normalize the gram matrix that will be returned
	
  unsigned int n = _collections.size();
	
  CImg<double> nGram (n,n);
	
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      if (i != j)
	nGram(i,j) = _gram(i,j)/sqrt(_gram(i,i)*_gram(j,j));
	
  for (unsigned int i=0; i<n; ++i)
    nGram (i,i)=1;
	
  return nGram;
}

void Dataset::showGramMatrix() const
{
  unsigned int n=_collections.size();
	
  for (unsigned int i=0; i<n; ++i)
    {
      cout << _parameters[i] << " 0:" << i+1;
		
      for (unsigned int j=0; j<n; ++j)
	{
	  cout.precision(DBL_DIG);
	  cout << " " << j+1 << ":" << _gram(i,j);
	}
		
      cout << endl;
    }
}

bool Dataset::isGramMatrixPD(){
  unsigned int n=_collections.size();
  assert (n > 0);
  CImg<double> eigvals;
  CImg<double> eigvects;
  _gram.symmetric_eigen(eigvals,eigvects);
  bool isPD = true;
  for (unsigned int i=0; i<n; ++i)
    {
      cout << eigvals[i] << endl;
      // cout << eigvects[i] << endl;
      //   if(eigvals[i] < 0)
      // 	return false;
      // 
      if(eigvals[i] < 0)
	isPD = false;
    }
  return isPD;
}

void Dataset::regularizeGramMatrix(){
  unsigned int n=_collections.size();
  assert (n > 0);
  CImg<double> eigvals;
  CImg<double> eigvects;
  _gram.symmetric_eigen(eigvals,eigvects);
  if(eigvals[n-1] < 0)
    for (unsigned int i=0; i<n; ++i)
      _gram(i,i) -= eigvals[n-1];
}


bool Dataset::isSymmetric(){
  unsigned int n=_collections.size();
  assert (n > 0);
  bool isSymmetric = true;
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      if(fabs(_gram(i,j) - _gram(j,i)) > DBL_EPSILON)
	isSymmetric = false;
  return isSymmetric;
}

void Dataset::showGramMatrixRaw() const
{
  unsigned int n=_collections.size();

   unsigned int toto;

  for (unsigned int i=0; i<n; ++i)
    {
       toto=0;
      for (unsigned int j=0; j<n; ++j)
	{
	  cout.precision(DBL_DIG);
	  cout <<  _gram(i,j)  << ",";
	   if(_gram(i,j)!=0)
	     toto++;
	}
      cout << endl;
       if(toto<2)
       	cout<<"Ligne "<<i<<" avec une valeur non nulle !"<<endl;
    }
}

void Dataset::showGramMatrixMatlab(char * filename) const
{
  ofstream outfile_gram (filename);

  unsigned int n=_collections.size();
  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int j=0; j<n; ++j)
	{
	  outfile_gram.precision(DBL_DIG);
	  outfile_gram <<  _gram(i,j);
	  outfile_gram.flush();
	  if(j!=n-1)
	    outfile_gram.flush() << ", ";
	}
      outfile_gram << endl;
    }
  outfile_gram.close();
}

void Dataset::delete_first ()
{
  unsigned int n = _collections.size();
	
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n-1; ++j)
      _gram(i, j)=_gram(i, j+1);
	
  for (unsigned int i=0; i<n-1; ++i)
    for (unsigned int j=0; j<n; ++j)
      _gram(i, j)=_gram(i+1, j);
	
  _collections.pop_front();
  _parameters.pop_front();
	
  _gram.resize(n, n);
}

void Dataset::delete_last ()
{
  unsigned int n = _collections.size();
	
  // 	for (unsigned int i=0; i<n; ++i)
  // 		for (unsigned int j=0; j<n-1; ++j)
  // 			_gram(i, j)=_gram(i, j+1);
	
  // 	for (unsigned int i=0; i<n-1; ++i)
  // 		for (unsigned int j=0; j<n; ++j)
  // 			_gram(i, j)=_gram(i+1, j);
	
  _collections.pop_back();
  _parameters.pop_back();
	
  _gram.resize(n, n);
}



void Dataset::add (Collection* m, double parameter, GraphKernel* kg)
{
  _collections.push_back(m);
  _parameters.push_back(parameter);
	
  if (kg != NULL)
    {
      unsigned int n = _collections.size();
		
      _gram.resize(n, n,-100,-100,0,0,0,0);

      for (unsigned int i=0; i<n; ++i)
	{
	  double x = (*kg) (_collections[n-1], _collections[i]);
			
	  _gram(n-1, i) = x;
	  _gram(i, n-1) = x;
	}
    }
}



int Dataset::find (pandore::Collection* col) const
{
  for (unsigned int i=0; i<_collections.size(); ++i)
    if (_collections[i] == col)
      return i;
	
  return -1;
}

void Dataset::destroy ()
{
  for (unsigned int i=0; i<_collections.size(); ++i)
    delete _collections[i];
}


// list<Dataset> Dataset::split(int nb_classes){
  
//   list<Dataset> l;
//   return l;
// }

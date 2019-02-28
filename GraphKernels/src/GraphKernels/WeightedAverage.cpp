/**
 * @file WeightedAverage.cpp
 *
 */

#include "WeightedAverage.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;
using namespace cimg_library;
using namespace pandore;

WeightedAverage::WeightedAverage (GraphKernel* kg, Dataset* dataset)
  {
    _kg = kg;
    _dataset = dataset;
  }


double WeightedAverage::operator() (Collection* col)
{
  unsigned int nDataset = _dataset->size();

  vector< pair<double,unsigned int > * > nonZ;

  for (unsigned int i=0; i<nDataset; ++i)
    {
      double k = (*_kg)(_dataset->getCollection(i), col);
      if(k!=0.0)
	nonZ.push_back(new pair<double,unsigned int>(k,i));
    }

  unsigned int n = nonZ.size();
  
  if(n==0)
    {
      cout<<"Cette molécule possède un stéreocentre unique !!"<<endl;
      return 0;
    }

  CImg<double> K (1,n);
  CImg<double> y (n,1);
	
  for (unsigned int j=0; j<n; ++j)
    {
      unsigned int i=nonZ[j]->second;
      double k = nonZ[j]->first;
      
      double a = (*_kg) (col, col);
      double b = (*_kg) (_dataset->getCollection(i),_dataset->getCollection(i));
      k /= sqrt(a*b);	
     
      K(0,j) = k;
      y(j,0) = _dataset->getParameter(i);
    }
	
  //CImg<double> gram = _dataset->getGramMatrix(true);


  //CImg<double> gram_nonZ (n,n);


  //CImg<double> tmp = gram_nonZ + CImg<double>(n,n).identity_matrix()*_lambda;
   

  //CImg<double> tmp2 = tmp.invert();
  
  //cout<<"N vaut:"<<n<<endl;
  CImg<double> res = y*K/(K.magnitude(1));//y*tmp2*K;


  /*if(n==2){
    cout<<"TESTTTTTT"<<endl;
    cout<<tmp(0,0)<<" "<<tmp(1,0)<<" "<<tmp(0,1)<<" "<<tmp(1,1)<<endl;
    cout<<"k1="<<K(0,0)<<" k2="<<K(0,1)<<endl;


    cout<<"y1="<<y(0,0)<<" y2="<<y(1,0)<<endl;
    cout<<"RES="<<res(0,0)<<endl;
    cout<<"lamda="<<_lambda<<endl;
    }*/

  return res(0,0);
}

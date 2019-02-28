/*
 * @file CombinedKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb 29 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "CombinedKernel.h"
#include <cmath>

using namespace std;
using namespace pandore;

CombinedKernel::CombinedKernel(vector<GraphKernel *> kernels,vector<double> weights,int typeOfCombination){
  _kernels = kernels;
  _weights = weights;
  _typeOfCombination = typeOfCombination;
}

double CombinedKernel::operator()(Collection* c1, Collection* c2){
  double k=0;
  if(_typeOfCombination==0)
    for(unsigned int i=0;i<_kernels.size();i++)
      k+=_weights[i]*(*_kernels[i])(c1,c2);
  if(_typeOfCombination==1)
    {
      k=1;
      for(unsigned int i=0;i<_kernels.size();i++)
	{
	  if((*_kernels[i])(c1,c2)!=0)
	    {
	      k*=pow((*_kernels[i])(c1,c2),_weights[i]);
	    }
	  else
	    {
	      k=0;
	    }
	}
    }
  return k;
}

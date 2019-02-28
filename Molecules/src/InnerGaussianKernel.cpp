/*
 * @file InnerGaussianKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <cmath>
#include <map>

#include "TreeletInnerGaussianKernel.h"

TreeletInnerGaussianKernel::TreeletInnerGaussianKernel(double sigma){
  _sigma = sigma;
}


double TreeletInnerGaussianKernel::operator() (treelet_spectrum spec_1, 
					       treelet_spectrum spec_2,
					       int treelet_type){
  double scalar = 0.0;
  
  treelet_spectrum::iterator it1 = spec_1.begin();
  while(it1 != spec_1.end())
    {
      if(isAdmitted(it1->first, treelet_type))
	if(spec_2.find(it1->first) != spec_2.end())
	  {
	    double diff = fabs(spec_2.find(it1->first)->second - it1->second);
	    double gaussian = exp(-(diff * diff)/(2*_sigma*_sigma));
	    double inner = spec_2.find(it1->first)->second * it1->second;
	    scalar += getWeight(it1->first,treelet_type)*inner*gaussian;
	  }
      it1++;
    }      
  return scalar;
}


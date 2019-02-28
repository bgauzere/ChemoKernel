/*
 * @file PolynomialKernel.cpp
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

#include "TreeletPolynomialKernel.h"

/*Noyau de Treelets calculant une somme de noyaux polynomiaux*/


TreeletPolynomialKernel::TreeletPolynomialKernel(int d){
  _d = d;
}


double TreeletPolynomialKernel::operator() (treelet_spectrum spec_1, 
				     treelet_spectrum spec_2,
				     int treelet_type){
  double scalar = 0.0;
  
  treelet_spectrum::iterator it1 = spec_1.begin();
  while(it1 != spec_1.end())
    {
      if(isAdmitted(it1->first, treelet_type))
	if(spec_2.find(it1->first) != spec_2.end())
	  {
	    scalar += getWeight(it1->first,treelet_type)*pow(spec_2.find(it1->first)->second * it1->second,_d);
	  }
      it1++;
    }      
  return scalar;
}


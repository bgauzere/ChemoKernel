/*
 * @file IntersectionKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Mar 15 2011
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

#include "TreeletIntersectionKernel.h"

using namespace std;
double TreeletIntersectionKernel::operator() (treelet_spectrum spec_1, 
					      treelet_spectrum spec_2,
					      int treelet_type){
  double scalar = 0.0;
  treelet_spectrum::iterator it1 = spec_1.begin();
  while(it1 != spec_1.end())
    {
      if(isAdmitted(it1->first, treelet_type))
	if(spec_2.find(it1->first) != spec_2.end())
	  {
	    double a = spec_2.find(it1->first)->second;
	    double b = it1->second;
	    scalar +=  getWeight(it1->first,treelet_type)*((a<b)?a:b);
	  }
      it1++;
    }      
  return scalar;
}


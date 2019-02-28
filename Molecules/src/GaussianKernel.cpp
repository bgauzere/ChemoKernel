/*
 * @file GaussianKernel.cpp
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

#include "GaussianKernel.h"

using namespace std;

GaussianKernel::GaussianKernel(double sigma){
  _sigma = sigma;
}

// // V2 : Compare spec_1 U spec_2
// double GaussianKernel::operator() (treelet_spectrum spec_1, 
// 				   treelet_spectrum spec_2){
//   double scalar = 0.0;
  
//   //Parcours de tous spec1, comptage treelets de spec1 seulement ET treelets communs
//   for(treelet_spectrum::iterator it1 = spec_1.begin();it1!= spec_1.end();it1++){
//     double u = it1->second;
//     double v = 0;//Chemin inexistant ds spec 2 par défaut
//     if(spec_2.find(it1->first) != spec_2.end()) //Chemin existant ds spec 2
//       v = spec_2.find(it1->first)->second;
//     double diff = fabs(u - v);
//     scalar += exp(-(diff * diff)/(2*_sigma*_sigma));
//   }


//   //Parcours des presents que dans spec2
//   for(treelet_spectrum::iterator it2 = spec_2.begin();it2!= spec_2.end();it2++){
//     if(spec_1.find(it2->first) == spec_1.end()) //Chemin inexistant dans spec 1
//       {
// 	double u = it2->second;
// 	double v = 0;//Chemin inexistant ds spec 1
// 	double diff = fabs(u - v);
// 	scalar += exp(-(diff * diff)/(2*_sigma*_sigma));
//       }
//   }
// }

// V1 : Compare spec_1 \intersection spec_2
double GaussianKernel::operator() (double x, double y){
  double diff = x-y;
  return exp(-(diff * diff)/(_sigma*_sigma));
}

/*
 * @file CompleteGaussianKernel.cpp
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
#include <cassert>

#include "TreeletCompleteGaussianKernel.h"

using namespace std;

TreeletCompleteGaussianKernel::TreeletCompleteGaussianKernel(double sigma){
  _sigma = sigma;
}

// V2 : Compare Tous les treelets du Dataset
double TreeletCompleteGaussianKernel::operator() (treelet_spectrum spec_1, 
					   treelet_spectrum spec_2, 
					   int treelet_type){
  double scalar = 0.0;
  for(unsigned int i = 0; i < _vit_list[treelet_type].size(); i++)
    {
      if(isAdmitted(_vit_list[treelet_type][i], treelet_type))
	{
	  double v = 0.0;//Chemin inexistant ds spec 2 par défaut
	  if(spec_2.find(_vit_list[treelet_type][i]) != spec_2.end()) //Chemin existant ds spec 2
	    {
	      v = spec_2.find(_vit_list[treelet_type][i])->second;
	    }
	  double u = 0.0;//Chemin inexistant ds spec 2 par défaut
	  if(spec_1.find(_vit_list[treelet_type][i]) != spec_1.end()) //Chemin existant ds spec 1
	    {
	      u = spec_1.find(_vit_list[treelet_type][i])->second;
	    }
	  double diff = fabs(u - v);
	  scalar += getWeight(_vit_list[treelet_type][i],treelet_type)*exp(-(diff * diff)/(2*_sigma*_sigma));
	}
    }
  return scalar;
}

// V2 : Compare spec_1 U spec_2
// double CompleteGaussianKernel::operator() (treelet_spectrum spec_1, 
// 					   treelet_spectrum spec_2, 
// 					   int treelet_type){
//   double scalar = 0.0;
//   treelet_spectrum::iterator it1 = spec_1.begin();
//   while(it1 != spec_1.end())
//     {
//       double u = it1->second;
//       double v = 0.0;//Chemin inexistant ds spec 2 par défaut
//       if(spec_2.find(it1->first) != spec_2.end()) //Chemin existant ds spec 2
// 	v = spec_2.find(it1->first)->second;
//       double diff = fabs(u - v);
//       scalar += exp(-(diff * diff)/(2*_sigma*_sigma));
//       it1++;
//     }  
  
//   it1 = spec_2.begin();
//   while(it1 != spec_2.end())
//     {
//       if(spec_1.find(it1->first) == spec_1.end())//Existant dans spec_2, inexistant dans spec 1
//   	  {
//   	    double u = it1->second;
//   	    double v = 0;//Chemin inexistant ds spec 2 par défaut
//   	    double diff = u;
//   	    scalar += exp(-(diff * diff)/(2*_sigma*_sigma));
//   	  }
//       it1++;
//     }

  

  
//   //qParcours de tous spec1, comptage treelets de spec1 seulement ET treelets communs
//   for(treelet_spectrum::iterator it1 = spec_1.begin();it1!= spec_1.end();it1++){
//     if(isAdmitted(it1->first,treelet_type))
//       {
//   	double u = it1->second;
//   	double v = 0;//Chemin inexistant ds spec 2 par défaut
//   	if(spec_2.find(it1->first) != spec_2.end()) //Chemin existant ds spec 2
//   	  v = spec_2.find(it1->first)->second;
//   	double diff = fabs(u - v);
//   	scalar += getWeight(it1->first,treelet_type)*exp(-(diff * diff)/(2*_sigma*_sigma));
//       }
//   }
  
  
//   //Parcours des presents que dans spec2
//   for(treelet_spectrum::iterator it2 = spec_2.begin();it2!= spec_2.end();it2++){
//     if(isAdmitted(it2->first,treelet_type))
//       {
//   	if(spec_1.find(it2->first) == spec_1.end()) //Chemin inexistant dans spec 1
//   	  {
//   	    double u = it2->second;
//   	    double v = 0;//Chemin inexistant ds spec 1
//   	    double diff = fabs(u - v);
//   	    scalar += getWeight(it2->first,treelet_type)*exp(-(diff * diff)/(2*_sigma*_sigma));
//   	  }
//       }
//   }
//   return scalar;
// }


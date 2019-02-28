/*
 * @file MinMaxKernel.cpp
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


#include <algorithm>
#include <cmath>

#include "MinMaxKernel.h"
#include "MoleculeGraph.h"
#include "InnerProductKernel.h"
#include "GaussianKernel.h"

using namespace std;
using namespace pandore;

double MinMaxKernel::operator() (Collection* c1, Collection* c2){
  treelet_spectrum ** spec_1 = (treelet_spectrum **) c1->GETVALUE("labeled_spectrum",Long);
  treelet_spectrum ** spec_2 = (treelet_spectrum **) c2->GETVALUE("labeled_spectrum",Long);
  double min_sum = 0.0;
  double max_sum = 0.0;
  
  for(int i = 0; i < 13; i++) //Parcours des 13 Treelets
    {
  
      //Parcours de tous spec1, comptage treelets de spec1 seulement ET treelets communs
      for(treelet_spectrum::iterator it1 = spec_1[i]->begin();it1!= spec_1[i]->end();it1++){
	double u = it1->second;
	double v = 0;//Chemin inexistant ds spec 2 par défaut
	if(spec_2[i]->find(it1->first) != spec_2[i]->end()) //Chemin existant ds spec 2
	  v = spec_2[i]->find(it1->first)->second;
	min_sum += (u<v)?u:v;
	max_sum += (u>v)?u:v;
      }

      //Parcours des presents que dans spec2
      for(treelet_spectrum::iterator it2 = spec_2[i]->begin();it2!= spec_2[i]->end();it2++){
	if(spec_1[i]->find(it2->first) == spec_1[i]->end()) //Chemin inexistant dans spec 1
	  {
	    double u = it2->second;
	    double v = 0;//Chemin inexistant ds spec 1
	    min_sum += (u<v)?u:v;
	    max_sum += (u>v)?u:v;
	  }
      }
    }
  return min_sum / max_sum;  
}

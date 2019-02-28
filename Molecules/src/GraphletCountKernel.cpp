/**
 * @file GraphFileKernel.cpp
 * XXX: A renommer
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.0.0 (2010-10-27)
 */

#include <algorithm>
#include <cmath>

#include "GraphletCountKernel.h"
#include "MoleculeGraph.h"

using namespace std;
using namespace cimg_library;
using namespace pandore;

GraphletCountKernel::GraphletCountKernel(double sigma, double mu){
  _sigma  = sigma;
  _mu  = mu;
  _graphlet_coeff = new double[SIZE_SPECTRUM];
  memset(_graphlet_coeff, 1, sizeof(double)*SIZE_SPECTRUM);
}

void GraphletCountKernel::setCoeffs(double * coeffs, int size){
  delete [] _graphlet_coeff;
  _graphlet_coeff = new double[size];
  for(int i = 0; i < size; i++)
    _graphlet_coeff[i] = pow(coeffs[i],_mu);
}


/*Gaussian Kernel*/
double GraphletCountKernel::operator() (Collection* c1, Collection* c2){
  double * m1 = c1->GETARRAY("spectrum", Double);
  // double nb_1 = c1->GETVALUE("nb_graphlets",Double);

  double * m2 = c2->GETARRAY("spectrum", Double);
  // double nb_2 = c2->GETVALUE("nb_graphlets",Double);
 
  unsigned int n = c1->GETARRAYSIZE("spectrum", Double);
  double scalar = 0;
  for(unsigned int i=0;i<13;i++){
    double diff = _graphlet_coeff[i] * (fabs(m1[i]-m2[i]));
    scalar += exp(-(diff * diff)/2*_sigma*_sigma);
  }

  return scalar;  
}

/*Polynomial Kernel*/
// double GraphletCountKernel::operator() (Collection* c1, Collection* c2){
//  double * m1 = c1->GETARRAY("spectrum", Double);
//  double nb_1 = c1->GETVALUE("nb_graphlets",Double);

//  double * m2 = c2->GETARRAY("spectrum", Double);
//  double nb_2 = c2->GETVALUE("nb_graphlets",Double);
 
//  unsigned int n = c1->GETARRAYSIZE("spectrum", Double);
//  double scalar = 0;
 
//  double c = 1.0;
//  //scalar += (nb_1*nb_2)/(fabs(m1[i]-m2[i])+1);
//  for(int i=0;i<n;i++){
//    scalar += pow(m1[i] * m2[i],_sigma) + c;
//  }
// return scalar;  
// }

/*Cosinus Kernel*/
// double GraphletCountKernel::operator() (Collection* c1, Collection* c2){
//  double * m1 = c1->GETARRAY("spectrum", Double);
//  double nb_1 = c1->GETVALUE("nb_graphlets",Double);

//  double * m2 = c2->GETARRAY("spectrum", Double);
//  double nb_2 = c2->GETVALUE("nb_graphlets",Double);
 
//  unsigned int n = c1->GETARRAYSIZE("spectrum", Double);
//  double scalar = 0;
//  //scalar += (nb_1*nb_2)/(fabs(m1[i]-m2[i])+1);
//  for(int i=0;i<n;i++){
//    scalar += m1[i] * m2[i];
//  }
 
//  return (scalar/(nb_1*nb_2));  
// }



/*Linear Kernel*/
// double GraphletCountKernel::operator() (Collection* c1, Collection* c2){
//  double * m1 = c1->GETARRAY("spectrum", Double);
//  double nb_1 = c1->GETVALUE("nb_graphlets",Double);

//  double * m2 = c2->GETARRAY("spectrum", Double);
//  double nb_2 = c2->GETVALUE("nb_graphlets",Double);
 
//  unsigned int n = c1->GETARRAYSIZE("spectrum", Double);
//  double scalar = 0;
//  //scalar += (nb_1*nb_2)/(fabs(m1[i]-m2[i])+1);
//  for(int i=0;i<n;i++){
//    double diff = fabs(m1[i]-m2[i]);
//    if(diff == 0)
//      scalar += 1;
//    else
//      scalar += 1/pow((diff*2),_sigma);
//  }
//  return scalar;  
// }


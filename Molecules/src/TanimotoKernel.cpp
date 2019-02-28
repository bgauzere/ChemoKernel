/*
 * @file TanimotoKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * @description : computes the Tanimoto Kernel (Ralaivola 2005) on the treelet spectrum
 * 
 */

#include <algorithm>
#include <cmath>

#include "TanimotoKernel.h"
#include "MoleculeGraph.h"
#include "TreeletInnerProductKernel.h"
#include "TreeletGaussianKernel.h"
#include "TreeletBinaryKernel.h"
#include "TreeletRandomKernel.h"

using namespace std;
using namespace pandore;

TanimotoKernel::TanimotoKernel(){
  k = new TreeletInnerProductKernel();
}


double TanimotoKernel::operator() (Collection* c1, Collection* c2){
  treelet_spectrum ** spec_1 = (treelet_spectrum **) c1->GETVALUE("labeled_spectrum",Long);
  treelet_spectrum ** spec_2 = (treelet_spectrum **) c2->GETVALUE("labeled_spectrum",Long);
  double k_uv = 0.0;
  double k_uu = 0.0;
  double k_vv = 0.0;
  // A verifier, pourquoi que 0 ?
  // for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
  //   {
  //    k_uv += (*k)(spec_1[i][0],spec_2[i][0]);
  //    k_uu += (*k)(spec_1[i][0],spec_1[i][0]);
  //    k_vv += (*k)(spec_2[i][0],spec_2[i][0]);    
  //   }
  // return k_uv / (k_uu + k_vv - k_uv);  
  return 0;
}

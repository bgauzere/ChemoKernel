/*
 * @file RandomKernel.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "TreeletRandomKernel.h"

using namespace std;
double TreeletRandomKernel::operator() (treelet_spectrum spec_1, 
					treelet_spectrum spec_2,
					int treelet_type){
  srand ( time(NULL) );
  return rand();
}


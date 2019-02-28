/*
 * @file InnerProductKernel.cpp
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

#include "InnerProductKernel.h"

double InnerProductKernel::operator() (double x, double y){
  return x*y;
}


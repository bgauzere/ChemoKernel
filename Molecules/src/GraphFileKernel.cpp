/**
 * @file GraphFileKernel.cpp
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.0.0
 */

#include <algorithm>
#include "GraphFileKernel.h"

// using fstream constructors.
#include <iostream>
#include <fstream>

using namespace std;
using namespace cimg_library;
using namespace pandore;

GraphFileKernel::GraphFileKernel(char * gram_file){
  ifstream filestr (gram_file, fstream::in);
  int w,h;
  filestr >> w >> h;
  Gram = CImg<double>(w,h,1,1,0.0);
  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
      {
	double tmp;
	filestr >> tmp;
	Gram(i,j) = tmp;
    }

  // for(int j=0; j<h; j++)
  //   {
  //     for(int i=0; i<w; i++)
  // 	cout << Gram(i,j)<<" ";
  //     cout << endl;
  //   }

  //for(int j=0; j<h; j++)
  //  cout << "(" << j <<","<< j << ") =" << Gram(j,j) << endl;  
}

double GraphFileKernel::operator()(pandore::Collection* c1, pandore::Collection* c2){
  long id_1 = c1->GETVALUE("id",Long);
  long id_2 = c2->GETVALUE("id",Long);
  return Gram(id_1,id_2);
}

void GraphFileKernel::normalize()
{
  unsigned int n = Gram.width();
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      if (i != j)
	Gram(i,j) /= sqrt(Gram(i,i)*Gram(j,j));
		
  for (unsigned int i=0; i<n; ++i)
    Gram (i,i)=1;
}

void GraphFileKernel::RBF(double delta)
{
  unsigned int n = Gram.width();

   for (unsigned int i=0; i<n; ++i)
	for (unsigned int j=0; j<n; ++j)
	  if (i != j)
	    Gram(i,j) = exp(- (Gram(i,i)+Gram(j,j)-2*Gram(i,j))/(2*delta*delta));
		
      for (unsigned int i=0; i<n; ++i)
	Gram (i,i)=1;
}

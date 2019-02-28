/*
 * @file VectorFileKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Jun 21 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <algorithm>
#include "VectorFileKernel.h"
// using fstream constructors.
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
using namespace cimg_library;
using namespace pandore;

VectorFileKernel::VectorFileKernel(char * vector_file, MoleculesDataset * dataset, int dimension, Kernel * k){
  this->dimension = dimension;
  this->dataset_size = dataset->size();
  ifstream filestr (vector_file, fstream::in);
  
  for(int i=0;i<dataset->size();i++)
    filenames.push_back(dataset->getMoleculeFilename(i));
    
  while (! filestr.eof()){
    string filename;
    filestr >> filename;
    if(! filename.empty()){
      
      double * vec = new double[dimension];
      for(int d=0;d<dimension;d++)
	filestr >> vec[d];
      vectors.insert(pair<string,double*>(filename.append(".ct"),vec));
      cout << filename << endl;
      // for(int i=0;i<dimension;i++)
      // 	cout << vec[i] << " ";
      // cout << endl;
    }    
  }
  this->k=k;
}


void VectorFileKernel::computeGramMatrices(char * gram_file){
  ofstream outfile_gram (gram_file);
  for(int m=0;m<dimension;m++)
    for(int i = 0;i < dataset_size;i++){
      map<string, double * >::iterator it1 = vectors.find(filenames[i]);
      double * vec_1 = it1->second;
      assert(it1 != vectors.end());
      for(int j = 0;j < dataset_size;j++){
	map<string, double * >::iterator it2 = vectors.find(filenames[j]);
	assert(it2 != vectors.end());
	double * vec_2 = it2->second;
	outfile_gram << (*k)(vec_1[m],vec_2[m]) << ",";
      }
      outfile_gram << endl;
    }
}

double VectorFileKernel::operator()(pandore::Collection* c1, pandore::Collection* c2){
  long id_1 = c1->GETVALUE("id",Long);
  long id_2 = c2->GETVALUE("id",Long);
  double sum=0.0;
  map<string, double * >::iterator it1 = vectors.find(filenames[id_1]);
  assert(it1 != vectors.end());
  
  map<string, double * >::iterator it2 = vectors.find(filenames[id_2]);
  assert(it2 != vectors.end());

  double * vec_1 = it1->second;
  double * vec_2 = it2->second;
  for(int i=0;i<dimension;i++)
    sum += (*k)(vec_1[i],vec_2[i]);
  return sum;
}


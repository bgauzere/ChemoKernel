/*
 * @file splitDataset.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Fri Feb 17 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <list>
#include <time.h>
#include "MoleculesDataset.h"

using namespace std;

void usage (char * s)
{
  cerr << "Usage : "<< s << "path/to/dataset/ dataset.ds nb_folds" << endl;
  cerr << "Split given dataset in nb_folds folds using random pick-up." << endl;
}

int main (int argc, char** argv)
{	
  
  int step = atoi(argv[3]);
  vector< vector<int> > index(step);
  set<int> used;
  MoleculesDataset * dataset = new MoleculesDataset(argv[1],argv[2]);
  int N = dataset->size();
  vector<int> nb_elmts(step);
  int reste = N - (step * (N/step));
  for(int i=0;i<step;i++){
    int nb_elmts = N/step;
    nb_elmts = (i<reste)?nb_elmts+1:nb_elmts;
    cout << nb_elmts << endl;
    for(int j=0;j<nb_elmts;j++){
      srand(time(NULL));
      int selected = rand()%N;
      while(used.find(selected) != used.end()){
	selected  = rand()%N;
      }
      used.insert(selected);
      index[i].push_back(selected);
    }
  }

  for(int k = 0; k < step; k++) {
    dataset = new MoleculesDataset(argv[1],argv[2]);
    list<int> test_set;
    list<string> test_set_filename;
    list<double> test_set_parameters;
    stringstream testset_filename;
    testset_filename << "testset_rand_" << k << ".ds";
    ofstream testset_file (testset_filename.str().c_str());

    for(int i = 0; i < index[k].size(); i++){
      int n = index[k][i];
      test_set.push_back(n);
      testset_file << dataset->getMoleculeFilename(n) << " "  << dataset->getParameter(n) << endl;
    }
    dataset->eraseSome(test_set);
	
    stringstream trainset_filename;
    trainset_filename << "trainset_rand_" << k << ".ds";
    ofstream trainset_file (trainset_filename.str().c_str());
    for(int i=0;i<dataset->size();i++){
      trainset_file << dataset->getMoleculeFilename(i) << " "  << dataset->getParameter(i) << endl;
    }
    
    trainset_file.close();
    testset_file.close();
  }
  return 0;
}

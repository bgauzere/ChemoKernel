/*
 * @file moleculesClassificationSlowAdd.cpp
 */

#include <iostream>

#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "GraphEditDistance.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"

#include "KEditDistance.h"

#include <SVM.h>

using namespace std;
using namespace pandore;

struct Options
{
  char* trainset_path;
  char* dataset_file;
  double sigma;
  double lambda;
  double regularization;
  double c;
  bool quietMode;
};

/**
 * Show the program's usage.
 */

void showUsage ()
{
  cerr << "Usage : moleculesClassification trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-c C : Set the SVM parameter C (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
}

/**
 * Read the command line's options.
 * 
 * @param argc The number of arguments.
 * @param argv The command line's arguments.
 * @param options The structure in which the parameters will be stored.
 */

void readOptions (int argc, char** argv, Options* options)
{
  if (argc < 3)
    {
      showUsage();
      exit(1);
    }
	
  // The default values are fixed
  options->sigma = 2.0;
  options->lambda = 1.0;
  options->regularization = 1;
  options->c = 1.0;
  options->quietMode = false;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];

  // Parsing the arguments
  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-l", 2) == 0)
	{
	  options->lambda = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else
	{
	  cerr << "Unknown parameter : " << argv[i] << endl;
	  showUsage();
	  exit(1);
	}
    }
}

int main (int argc, char** argv)
{	
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  Options options;
  readOptions(argc, argv, &options);
	
  GraphEditDistance* edit = new MoleculeGraphEditDistance;
  MoleculeGraph::initTable();
  double accuracy = 0.0;
	
  unsigned int N = 68;
	
  for (unsigned int i=0; i<N; ++i)
    {
      MoleculesDataset trainset (options.trainset_path, options.dataset_file);
      double classe = trainset.getParameter(i);
      Collection * to_test  = trainset.getCollection(i);
      list<int> to_erase;
      to_erase.push_back(i);
      trainset.eraseSome(to_erase);
      
      GraphKernel*  kgraph = new LaplacianKernel (edit, trainset, options.sigma, options.regularization, options.lambda);
      trainset.computeGramMatrix(kgraph, true);
      
      SVM svm (trainset, kgraph, SVM::C_SVC, options.c);
      
      double c = classe;
      double x = svm.predict(to_test);
      if(!options.quietMode)
	cout << c << " --> " << x << endl;
      
      if (c == x)
	++accuracy;
      
      trainset.destroy();
      delete kgraph;
    }
  
  if(!options.quietMode)
    cout << endl;
  
  cout << "Accuracy : " << 100*accuracy/N << "% (" << accuracy << " / " << N << ")" << endl;

  delete edit;
	
  return 0;
}

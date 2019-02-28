/*
 * @file moleculesRegressionDirectory.cpp
 *
 * @description This program performs a Kernel Ridge regression on gram matrix files present in a directory.
 *
 * @author Benoit Gauzere <benoit.gauzere@ensicaen.fr>
 */

#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <pandore.h>

#include "CImg.h"
#include "MoleculesDataset.h"
#include "MoleculeGraph.h"
#include "GraphFileKernel.h"

#include "KernelRidge.h"

using namespace std;
using namespace pandore;

void leaveOneOut (MoleculesDataset & dataset, GraphFileKernel* kgraph, double alpha, bool quietMode)
{	
  KRegression* kr = new KernelRidge (kgraph, &dataset, alpha);
  int N = dataset.size();

  //Stats data
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;
  double * bp_exps = new double[N];
  for (int i=0; i<N; ++i)
    {
      double bp_exp = dataset.getParameter(0);
      Collection* col = dataset.getCollection(0);
      dataset.delete_first();
      double bp = (*kr) (col);
      double err = bp_exp - bp;
      err_quad += err * err;
      if(!quietMode){
	cout << "Molecule " << i << " : " << bp << " (" << bp_exp << ")";
	cout << "\t --> \t" << abs(err)  << "°" << endl;
      }
      err_moy += abs(bp_exp - bp);
      bp_mean += bp_exp; 
      bp_exps[i] = bp_exp;		
      dataset.add (col, bp_exp, kgraph);
    }

  bp_mean /= N;
  std_dev = sqrt(err_quad / N); 
  err_moy /= N;

  double r_den = 0.0;
  for(int i =0;i < N; i++){
    r_den += pow(bp_exps[i] - bp_mean,2);
  }
  double r = sqrt(1- (err_quad) / (r_den));
  cout << "\t Average error : " << err_moy  << endl;
  cout << "\t Standard Deviation : " << std_dev << endl;
  cout << "\t R : " << r << endl;
  delete [] bp_exps;
  delete kr;
}

struct Options
{
  char* trainset_path;
  char* dataset_file;
  char * directory_gram;
  double alpha;
  bool quietMode;
};

void usage(char * s)
{
  cerr << "Usage : " << s << " trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-k : path to gram matrices directory" << endl;
  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
}

void readOptions (int argc, char** argv, Options* options)
{
  cout << endl;
  if (argc < 3)
    {
      usage(argv[0]);
      exit(1);
    }
	
  options->directory_gram = "";
  options->alpha = 1.0;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->quietMode = false;

  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->directory_gram = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-a", 2) == 0)
	{
	  options->alpha = atof(argv[i+1]);
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
	  usage(argv[0]);
	  exit(1);
	}
    }
}

int main (int argc, char** argv)
{	
  //Print command for archiving results
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }

  Options options;
  readOptions(argc, argv, &options);
  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.trainset_path, options.dataset_file);
  if(!options.quietMode)
    {
      cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
      cout << "Number of molecules : " << dataset.size() << endl;
      cout << "Alpha : " << options.alpha << endl;
      cout << endl;
    }
  DIR * dir = opendir(options.directory_gram);
  struct dirent * next_entry;
  while((next_entry = readdir(dir)) != NULL){
    //Regression on each file of dir
    if(next_entry->d_name[0] != '.'){
      char * filename = new char[strlen(options.directory_gram)+strlen(next_entry->d_name)+1];
      filename = strcpy(filename, options.directory_gram);
      filename = strcat(filename, "/");
      filename = strcat(filename,next_entry->d_name);
      cout << "--------------------------------" << endl;
      cout << next_entry->d_name << " : " << endl;
      GraphFileKernel * kgraph = new GraphFileKernel(filename);
      dataset.computeGramMatrix(kgraph,false);
      leaveOneOut(dataset, kgraph, options.alpha,options.quietMode);
      cout << "--------------------------------" << endl;
      delete kgraph;
      delete [] filename;
    }
  }

  dataset.destroy();
  return 0;
}

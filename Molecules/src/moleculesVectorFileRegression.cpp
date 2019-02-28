/*
 * @file moleculesVectorFileRegression.cpp
 *
 * @description 
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
#include "VectorFileKernel.h"
#include "IntersectionKernel.h"
#include "GaussianKernel.h"
#include "InnerProductKernel.h"

using namespace std;
using namespace pandore;

double predict (MoleculesDataset & trainset, GraphKernel* kgraph, double lambda, Collection* col)
{	
  KRegression* kr = new KernelRidge (kgraph, &trainset, lambda);
  return (*kr) (col);
}


void leaveOneOut (MoleculesDataset & dataset, GraphKernel* kgraph, double alpha, bool quietMode)
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
      char * filename = dataset.getMoleculeFilename(i);

      dataset.delete_first();
      double bp = (*kr) (col);
      double err = bp_exp - bp;
      err_quad += err * err;
      if(!quietMode){
	cout << "Molecule " << i << "(" << filename << ") : " << bp << " (" << bp_exp << ")";
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
  char* dataset_path;
  char* dataset_file;
  char * vector_file;
  char * gram_file;
  double alpha;
  int useTestSet;
  int dim;
  int K;
  double sigma;
};

void usage(char * s)
{
  cerr << "Usage : " << s << " dataset_path dataset_file" << endl;
  cerr << "options:" << endl;
  cerr << "-v : path to vector file" << endl;
  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-t nb : Use 1 molecule / nb as test, others as training set" << endl;
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Polynomial Kernel" << endl;
  cerr << "-d dim : Vectors dimension" << endl;
  cerr << "-s sigma : Set the parameter for gaussian kernel or power for polynomial kernel" << endl;
  cerr << "-m filename : write Gram matrix in Matlab format into filename" << endl;

}

void readOptions (int argc, char** argv, Options* options)
{
  cout << endl;
  if (argc < 3)
    {
      usage(argv[0]);
      exit(1);
    }
	
  options->vector_file = "";
  options->alpha = 1.0;
  options->dataset_path = argv[1];
  options->dataset_file = argv[2];
  options->useTestSet = 0;
  options->sigma = 2.0;
  options->dim =1;
  options->gram_file = "";

  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-v", 2) == 0)
	{
	  options->vector_file = argv[i+1];
	  i+=2;
	}
      if (strncmp(argv[i], "-m", 2) == 0)
	{
	  options->gram_file = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-a", 2) == 0)
	{
	  options->alpha = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->useTestSet = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->K = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->dim = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
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
  cout << endl;
  Options options;
  readOptions(argc, argv, &options);
  MoleculeGraph::initTable();
  MoleculesDataset * dataset = new MoleculesDataset(options.dataset_path, options.dataset_file);
  
  Kernel *k=NULL;
  switch(options.K){
  case Kernel::IntersectionKernelType:
    k=new IntersectionKernel();
    break;
  case Kernel::GaussianKernelType:
    k=new GaussianKernel(options.sigma);
    break;
  case Kernel::InnerProductKernelType:
    k=new InnerProductKernel();
    break;
  default:
    cerr << "Kernel undefined" << endl;
    exit(EXIT_FAILURE);
  }

  
  VectorFileKernel* kgraph =  new VectorFileKernel(options.vector_file,dataset,options.dim,k);
  dataset->computeGramMatrix(kgraph,false);
  if(strlen(options.gram_file) > 0)
    kgraph->computeGramMatrices(options.gram_file);
  if(options.useTestSet == 0)
    leaveOneOut(*dataset, kgraph, options.alpha,0);
  else{ 
    //Stats data
    int N = dataset->size();
    double err_moy = 0;
    double std_dev = 0.0;
    double err_quad = 0.0;
    double * err_quad_by_struct = new double[SIZE_SPECTRUM];
    double bp_mean = 0.0;
    double * bp_exps = new double[N];
    int nb_molecules = 0;
    for(int i = 0; i < N; i ++){
      bp_exps[i] = dataset->getParameter(i);
    }
    int step = options.useTestSet;
    for(int k = 0; k < step; k++) {
      cout << k << "/" << step<< endl;
      //Dataset split
      list<int> to_test;
      list<Collection *> to_test_Collections;
      list<double> to_test_Parameters;
      for(int n = k; n < N; n=n+step){
	to_test.push_back(n);
	to_test_Collections.push_back(dataset->getCollection(n));
	to_test_Parameters.push_back(dataset->getParameter(n));
      }
      dataset->eraseSome(to_test);
      //The Dataset is splitted
      double * coeffs;
      dataset->computeGramMatrix(kgraph, false);
      
      //Test molecules prediction
      list<int>::iterator it = to_test.begin();
      list<Collection*>::iterator it_col = to_test_Collections.begin();    
      list<double>::iterator it_param = to_test_Parameters.begin();
      
      for(;it != to_test.end();++it, ++it_col, ++it_param){
	double bp = predict(*dataset, kgraph, options.alpha,*it_col);
	double bp_exp = *it_param;
	double err =  bp_exp - bp;
	err_quad +=  err * err;
	err_moy += abs(bp_exp - bp);
	bp_mean += bp_exp; 
	nb_molecules ++;
      }
      //Dataset cleaning
      delete dataset;
      dataset = new MoleculesDataset (options.dataset_path, options.dataset_file);

    }
    bp_mean /= nb_molecules;
    std_dev = sqrt(err_quad / nb_molecules); 
    err_moy /= nb_molecules;

    double r_den = 0.0;
    for(int i =0;i < N; i++){
      r_den += pow(bp_exps[i] - bp_mean,2);
    }
    double r = sqrt(1- (err_quad) / (r_den));

    //Print Results
    cout << endl;
    cout << "Average error : " << err_moy  << endl;
    cout << "Standard Deviation : " << std_dev << endl;
    cout << "R : " << r << endl;
  }
  
  delete kgraph;
  dataset->destroy();
  return 0;
}

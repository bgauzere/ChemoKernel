/*
 * @file moleculesRegressionCluster.cpp
 *
 * This program performs a Kernel Ridge regression on molecules datasets divided in clusters.
 *
 * Usage : moleculesRegressionCluster trainset_path dataset_file [options]
 * options:
 * -k kernel_type : Set the type of kernel function (default 0)
 *		0 -- Kernel Mean
 *		1 -- Kernel Weighted Mean
 *		2 -- Graph Laplacian Kernel
 *		3 -- Gaussian Kernel
 *              4 -- GraphletCountKernel
 * -a alpha : Set the parameter of the regression (default 1.0)
 * -s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel (default 2.0)
 * -r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)
 * -l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)
 * -f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)
 * -m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)
 * -e editdistance : Set the Edit Distance computing : 
 *              0 -- V1
 *              1 -- V2
 * -q : Quiet mode, display only the percentage of good classification
 * -n mu : Parametre mu pour le GraphletCountKernel
 * -c cluster_mode : Sets the way used to clustering the dataset
 *              0 -- By size
 *              
 * 
 * @author Benoit Gauzere <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-10)
 */

#include <iostream>
#include <vector>

#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "TrailKernel.h"
#include "ShapeBagTrails.h" 
#include <pandore.h>

#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "KRegression.h"
#include "GraphletCountKernel.h"

#include "KernelRidge.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"
#include "KEditDistance.h"

using namespace std;
using namespace pandore;

void leaveOneOut (MoleculesDataset & dataset, GraphKernel* kgraph, double lambda,bool quietMode)
{	
  KRegression* kr = new KernelRidge (kgraph, &dataset, lambda);
	
  double err_moy = 0;
	
  int N = dataset.size();
	
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
		
      //trop dependant de la temp d ebullition (grosse erreur sur temp proche de 0)
      // cout << "Molecule " << i << " : " << bp << " (" << bp_exp << ")";
      // cout << " --> " << abs(bp_exp - bp)/bp_exp*100 << "% (lambda = " << lambda << ")" << endl;
      // err_moy += (bp_exp == 0)?bp:abs(abs(bp_exp - bp)/bp_exp*100);
		
      double err = bp_exp - bp;
      err_quad += err * err;
		
      if(!quietMode)
	{
	  cout << "Molecule " << i << " : " << bp << " (" << bp_exp << ")";
	  cout << "\t --> \t" << abs(err)  << "° (lambda = " << lambda << ")" << endl;
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


  cout << endl;
  cout << "Average error : " << err_moy  << endl;
  cout << "Standard Deviation : " << std_dev << endl;
  cout << "R : " << r << endl;

  delete [] bp_exps;
  delete kr;
}

struct Options
{
  char* trainset_path;
  char* dataset_file;
  int kernel;
  double sigma;
  int regularization;
  double lambda;
  double alpha;
  double keptTrails;
  double trailsLength;
  int editdistance;
  bool quietMode;
  int clusterMode;
  double mu;

};

/**
 * Show the program's usage.
 */

void showUsage ()
{
  cerr << "Usage : moleculesRegression trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
	
  cerr << "-k kernel_type : Set the type of kernel function (default 0)" << endl;
  cerr << "\t0 -- Kernel Mean" << endl;
  cerr << "\t1 -- Kernel Weighted Mean" << endl;
  cerr << "\t2 -- Graph Laplacian Kernel" << endl;
  cerr << "\t3 -- Gaussian Kernel" << endl;
  cerr << "\t4 -- Graphlet Count Kernel" << endl;
	
  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel & Graphlet Count Kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-e editdistance : Set the Edit Distance computing : " << endl;
  cerr << "\t0 -- V1" << endl;
  cerr << "\t1 -- V2" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-c cluster_mode : Sets the way used to clustering the dataset" << endl;
  cerr << "\t0 -- By Size" << endl;
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
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  if (argc < 3)
    {
      showUsage();
      exit(1);
    }
	
  // The default values are fixed
	
  options->kernel = 0;
  options->sigma = 2.0;
  options->lambda = 1.0;
  options->regularization = 0;
  options->alpha = 1.0;
  options->keptTrails = 1.0;
  options->trailsLength = 15;
  options->editdistance = 0;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->quietMode = false;
  options->clusterMode = 0;

  // Parsing the arguments
	
  int i=3;
	
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->kernel = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-r", 2) == 0)
	{
	  options->regularization = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-l", 2) == 0)
	{
	  options->lambda = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-a", 2) == 0)
	{
	  options->alpha = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-f", 2) == 0)
	{
	  options->keptTrails = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-m", 2) == 0)
	{
	  options->trailsLength = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-e", 2) == 0)
	{
	  options->editdistance = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-n", 2) == 0)
	{
	  options->mu = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->clusterMode = atof(argv[i+1]);
	  i+=2;
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
  Options options;
  readOptions(argc, argv, &options);
	
  ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
	
  GraphEditDistance* edit = NULL;
  if(options.editdistance == 0){
    edit = new MoleculeGraphEditDistance;
  }else{
    //V2
    edit = new MoleculeGraphEditDistanceV2;
  }
	
  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.trainset_path, options.dataset_file);
  vector<MoleculesDataset*> datasets = dataset.splitBySize();
  vector<GraphKernel*> kgraphs(datasets.size());
  double * coeffs;
	
  switch(options.kernel)
    {
    case 0:
      cout << "Kernel : KMean" << endl;
      for(int i =0;i<datasets.size();i++)
	kgraphs[i] = new KMean(sbt, NULL);
      break;
			
    case 1:
      cout << "Kernel : KWMean" << endl;
      cout << "Sigma : " << options.sigma << endl;
      for(int i =0;i<datasets.size();i++)
	kgraphs[i] = new KWMean (sbt, NULL, options.sigma);
      break;
			
    case 2:
      cout << "Kernel : Graph Laplacian Kernel" << endl;
      cout << "Regularization type : ";
			
      switch (options.regularization)
	{
	case 0: cout << "No regularization" << endl; break;
	case 1: cout << "Regularized Laplacian" << endl; break;
	case 2: cout << "Diffusion Process" << endl; break;
	case 3: cout << "Inverse Cosine" << endl; break;
	case 4: cout << "One-Step Random Walk" << endl; break;
	}
			
      cout << "Sigma : " << options.sigma << endl;
      cout << "Lambda : " << options.lambda << endl;
      
      for(int i =0;i<datasets.size();i++)
	kgraphs[i] = new LaplacianKernel (edit, *datasets[i], options.sigma, options.regularization, options.lambda);
      cout << "Poids :" << endl;
      //((LaplacianKernel *)kgraph)->printWeightMatrix();
      break;
			
    case 3:
      for(int i =0;i<datasets.size();i++)
	kgraphs[i] = new KEditDistance (edit, options.sigma);
      break;

    case 4:  
      for(int i =0;i<datasets.size();i++)
	{
	  datasets[i]->computeSpectrums();
	  datasets[i]->computeGraphletCorrelation();
	  kgraphs[i] = new GraphletCountKernel(options.sigma,options.mu);
	  
	  coeffs = new double[SIZE_SPECTRUM];
	  for(int k = 0;k < SIZE_SPECTRUM;k++)
	    coeffs[k] = datasets[i]->getCorrelationCoeff(k);
	  ((GraphletCountKernel*)kgraphs[i])->setCoeffs(coeffs,SIZE_SPECTRUM);
	}
      break;
    default:
      cerr << "Error : Invalid kernel" << endl;
      exit(1);
    }

  for(int i=0;i<datasets.size();i++){
    if(datasets[i]->size() > 5)
      {
	

	if(!options.quietMode){
	  cout << "Dataset " << i << " : " << options.trainset_path << "/" << options.dataset_file << endl;
	  cout << "Number of molecules : " << datasets[i]->size() << endl;
	  cout << "Alpha : " << options.alpha << endl;
	  cout << endl;
	}
	
	// We compute the gram matrix of the training set
	
	datasets[i]->computeGramMatrix(kgraphs[i], false);
	
	// Start the leave-one-out
	
	leaveOneOut(*datasets[i], kgraphs[i], options.alpha,options.quietMode);
	
	// Destruction of the training set
	
	datasets[i]->destroy();
	delete kgraphs[i];
      }
  }
	
  delete edit;
	
  return 0;
}

/*
 * 
 *
 * This program performs a Kernel Ridge regression on molecules datasets.
 * Version moins Naive : on calcule le dataset une fois pour toutes
 *
 * Usage : moleculesRegression trainset_path dataset_file [options]
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
 * 
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
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

void leaveOneOut (MoleculesDataset & dataset, double lambda,bool quietMode,
		  double sigma,   GraphEditDistance* edit)
{	
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;

  vector<Collection *> to_test;
  vector<double> bps_test;
  list<int> to_erase;

  //on enleve 10 % du dataset
  for(double i = 0; i < dataset.size();i+=(double)dataset.size()/10.){
    bps_test.push_back(dataset.getParameter(i));
    to_test.push_back(dataset.getCollection(i));
    to_erase.push_back(i);
  }
  dataset.eraseSome(to_erase);
  LaplacianKernel* kgraph = new LaplacianKernel (edit, dataset, sigma, 1.0, 1.0);
  dataset.computeGramMatrix(kgraph, false);
  KRegression* kr = new KernelRidge (kgraph, &dataset, lambda);

  
  int N = to_test.size();
  cout << N << endl;
  double * bp_exps = new double[N];
  
  for (int i=0; i<N; ++i)
    {
      cout << i << endl;
      double bp_exp = bps_test[i];
      Collection* col = to_test[i];
      //kgraph_current->fastAdd(col); //K est calculé
      //XXX : Dataset a mettre a jour
      //dataset.add(col,bp_exp,kgraph_current); //MoleculeGraph non utilisé dans KernelRidge
      //dataset.computeGramMatrix(kgraph_current);//On copie K ds gram
      
      double bp = (*kr) (col);
		
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
      // delete kgraph_current;
      // delete kr;
      // dataset.delete_last();
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
}


struct Options
{
  char* trainset_path;
  char* dataset_file;
  double sigma;
  double alpha;
  int editdistance;
  bool quietMode;
};

void showUsage ()
{
  cerr << "Usage : moleculesRegression trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
	
  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel & Graphlet Count Kernel (default 2.0)" << endl;
  cerr << "-e editdistance : Set the Edit Distance computing : " << endl;
  cerr << "\t0 -- V1" << endl;
  cerr << "\t1 -- V2" << endl;
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
	
  options->sigma = 2.0;
  options->alpha = 1.0;
  options->editdistance = 0;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->quietMode = false;

  // Parsing the arguments
	
  int i=3;
	
  while (i<argc)
    {
      if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-a", 2) == 0)
	{
	  options->alpha = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-e", 2) == 0)
	{
	  options->editdistance = atoi(argv[i+1]);
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
  Options options;
  readOptions(argc, argv, &options);
	
  
  GraphEditDistance* edit = NULL;
  if(options.editdistance == 0){
    edit = new MoleculeGraphEditDistance;
  }else{
    //V2
    edit = new MoleculeGraphEditDistanceV2;
  }
	
  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.trainset_path, options.dataset_file);
	
  cout << "Kernel : Regularized Graph Laplacian Kernel" << endl;
  cout << "Sigma : " << options.sigma << endl;
  

  if(!options.quietMode){
    cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
    cout << "Number of molecules : " << dataset.size() << endl;
    cout << "Alpha : " << options.alpha << endl;
    cout << endl;
  }
	
  leaveOneOut(dataset, options.alpha,options.quietMode, options.sigma,edit);
  dataset.destroy();
  delete edit;
  	
  return 0;
}

/*
 * @file moleculesClassificationKPCA.cpp
 *
 * @version 1.1.0 (2010-07-10)
 */

#include <iostream>
#include <cfloat>
#include "CImg.h"

#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "TrailKernel.h"
#include "ShapeBagTrails.h"
#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "GraphEditDistance.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"
#include "GraphletCountKernel.h"
#include "Kashima.h"
#include "KEditDistance.h"
#include "TreeletKernel.h"
#include "KernelPCA.h"
#include "TreeletType.h"

#define DIMENSION 3

using namespace std;
using namespace pandore;
using namespace cimg_library;
struct Options
{
  char* trainset_path;
  char* dataset_file;
  char* testset_path;
  char* testset_dataset_file;
  int kernel;
  double sigma;
  int regularization;
  double lambda;
  double keptTrails;
  double trailsLength;
  bool useTestset;
  bool quietMode;
  int customSizeTrainset;
  int dimension;
  bool normalize;
  KernelType spectrumKernel;
};

/**
 * Show the program's usage.
 */

void showUsage ()
{
  cerr << "Usage : moleculesClassification trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
	
  cerr << "-k kernel_type : Set the type of kernel function (default 0)" << endl;
  cerr << "\t0 -- Kernel Mean" << endl;
  cerr << "\t1 -- Kernel Weighted Mean" << endl;
  cerr << "\t2 -- Graph Laplacian Kernel" << endl;
  cerr << "\t3 -- Gaussian Kernel" << endl;
  cerr << "\t4 -- Graphlet Count Kernel" << endl;
  cerr << "\t5 -- Random Walks Kernel" << endl;
  cerr << "\t6 -- Treelet Kernel" << endl;
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;

  cerr << "-t testset_path dataset_file : Set the test set to use (default none)" << endl;
  cerr << "-c C : Set the SVM parameter C (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-n : Size of the train dataset considered to predict the class of the molecule tested" << endl;
  cerr << "-N : Normalize Treelet Kernel" << endl;
  cerr << "-d : Dimension of the PCA" << endl;
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
	
  options->kernel = 0;
  options->sigma = 2.0;
  options->lambda = 1.0;
  options->regularization = 0;
  options->keptTrails = 1.0;
  options->trailsLength = 15;
  options->useTestset = false;
  options->quietMode = false;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->customSizeTrainset = 0;
  options->dimension = DIMENSION;
  options->normalize = true;
  options->spectrumKernel = InnerProductKernelType;
  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->kernel = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->testset_path = argv[i+1];
	  options->testset_dataset_file = argv[i+2];
	  options->useTestset = true;
	  i+=3;
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
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-n", 2) == 0)
	{
	  options->customSizeTrainset = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->dimension = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-N", 2) == 0)
	{
	  options->normalize = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->spectrumKernel = (KernelType)(atoi(argv[i+1]));
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
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  Options options;
  readOptions(argc, argv, &options);
	
  ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
	
  GraphKernel* kgraph = NULL;

  GraphEditDistance* edit = new MoleculeGraphEditDistance;
	
  // Creation of the training set

  MoleculeGraph::initTable();
  //trainset : Base d'apprentissage
  MoleculesDataset trainset (options.trainset_path, options.dataset_file);
  if(options.customSizeTrainset == 0) 
    options.customSizeTrainset = trainset.size();

  
  //testset : Base de test
  MoleculesDataset testset;
	
  if (options.useTestset){
    testset.loadDataset (options.testset_path, options.testset_dataset_file);
    //Reduction de la taille de trainset
    //   while(trainset.size()>options.customSizeTrainset)
    //       trainset.delete_first();
    trainset.reduceToN(options.customSizeTrainset);
    cout << "trainset Reduce " << trainset.size() << endl; 
  }
  
  // Kernels initializations
	
  switch(options.kernel)
    {
    case 0:
      if(!options. quietMode)
	cout << "Kernel : KMean" << endl;
      kgraph = new KMean(sbt, new Kashima());
      break;
			
    case 1:
      if(!options.quietMode){
	cout << "Kernel : KWMean" << endl;
	cout << "Sigma : " << options.sigma << endl;
      }
			
      kgraph = new KWMean (sbt, NULL, options.sigma);
      break;
			
    case 2:
      if(!options.quietMode){
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
      }	

      if (options.useTestset)
	kgraph = new LaplacianKernel (edit, trainset, testset, options.sigma, options.regularization, options.lambda);
      else
	kgraph = new LaplacianKernel (edit, trainset, options.sigma, options.regularization, options.lambda);
      break;
			
    case 3:
      	cout << "Kernel : Edit Distance" << endl;
      kgraph = new KEditDistance (edit, options.sigma);
      break;	
    case 4:
      cout << "Kernel : Unlabeled Treelet Kernel" << endl;
      trainset.computeSpectrums();
      trainset.printSpectrums();
      kgraph = new GraphletCountKernel();
      break;
    case 5:
      cout << "Kernel : Random Walks Kernel" << endl;
      kgraph = new RandomWalkKernel(new Kashima(),1.0,1.0);
      break;			
    case 6:  
      cout << "Kernel : Treelet Kernel" << endl;      
      cout << "Spectrum Kernel : ";
      switch(options.spectrumKernel)
	{
	case IntersectionKernelType:
	  cout << "Intersection Kernel";
	  break;
	case GaussianKernelType:
	  cout << "Gaussian Kernel";
	  break;
	case InnerProductKernelType:
	  cout << "Inner Product Kernel";
	  break;
	case BinaryKernelType:
	  cout << "Binary Kernel";
	  break;
	case RandomKernelType:
	  cout << "Random Kernel";
	  break;
	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      trainset.computeLabeledSpectrums();
      if(options.normalize)
	trainset.normalizeLabeledSpectrums();
      kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
      break;
    default:
      cerr << "Error : Invalid kernel" << endl;
      exit(1);
    }
	
  if(!options.quietMode)
    {
      cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
      cout << "Number of molecules : " << trainset.size() << endl;
      cout << endl;
	
    }

  // We compute the Gram matrix of the training set
  trainset.computeGramMatrix(kgraph, true);
  KernelPCA * pca_projection = new KernelPCA(kgraph, trainset);
  
  double ** projection_result = pca_projection->projectTrainingSet(options.dimension);  

  //on cherche le plus proche graphe en distance euclidienne  (2D)
  double accuracy = 0.0;
  unsigned int N = trainset.size();
  
  for (unsigned int i=0; i<N; ++i)
    {
      int j_min = -1;
      double distance_min = DBL_MAX;
      for(unsigned int j=0; j<N; ++j)
	{
	  if(i!=j)
	    {
	      double distance = 0.0;
	      for(int dimension = 0 ; dimension < options.dimension; dimension++)
		distance += pow(projection_result[i][dimension] - projection_result[j][dimension],2);
	      
	      if(distance < distance_min)
		{
		  distance_min = distance;
		  j_min = j;
		}
	    }
	}
      
      double x = trainset.getParameter(j_min);      
      if(!options.quietMode)
	{
	  cout << trainset.getParameter(i) << " --> " << x << endl;
	}
      
      if (trainset.getParameter(i) == x)
	++accuracy;
    }
      

  if(!options.quietMode)
    cout << endl;

  if(options.quietMode)
    cout << options.customSizeTrainset << " " << 100*accuracy/N << endl;
  else
    cout << "Accuracy : " << 100*accuracy/N << "% (" << accuracy << " / " << N << ")" << endl;

  trainset.destroy();
	
  delete edit;
  delete kgraph;
	
  return 0;
}

/*
 * @file moleculesKPCA.cpp
 *
 * @description This program computes a formatted dlm file (format: x y value) 
 *                   or plots the PCA of the gram matrix (2D). 
 *              The file should be used with cplot G'MIC's function
 * 
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-10)
 *   TODO : Rajouter autres kernels
 */

#include <iostream>
#include <fstream>
#include <cfloat>
#include <pandore.h>

//Kernel includes
#include "GraphKernel.h"
#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "GraphletCountKernel.h"
#include "LaplacianKernel.h"
#include "KEditDistance.h"
#include "TrailKernel.h"
#include "ShapeBagTrails.h" 
#include "TreeletKernel.h"
#include "string_utils.h"
//#include "Kashima.h"
#include "TreeletType.h"

#include "MoleculesDataset.h"
#include "MoleculeGraph.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"

#include "KernelPCA.h"

//For instant plotting
extern "C" {
#include "gnuplot_i.h"
}


using namespace std;
using namespace pandore;

struct Options
{
  char* dataset_path;
  char* dataset_file;
  int kernel;

  double keptTrails;
  double trailsLength;
  double sigma;
  int regularization;
  double lambda;
  double alpha;
  int editdistance;
  double mu;
  bool quietMode;
  bool normalize;
  KernelType spectrumKernel;
  char* output_file;
};

void usage (char * s)
{
  cerr << "Usage : " << s << " dataset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;	
  cerr << "-k kernel_type : Set the type of kernel function (default 0)" << endl;
  cerr << "\t0 -- Kernel Mean" << endl;
  cerr << "\t1 -- Kernel Weighted Mean" << endl;
  cerr << "\t2 -- Graph Laplacian Kernel" << endl;
  cerr << "\t3 -- Gaussian Kernel" << endl;
  cerr << "\t4 -- Unlabeled Treelet Kernel" << endl;
  cerr << "\t5 -- Randow Walks Kernel" << endl;  
  cerr << "\t6 -- Treelet Kernel" << endl;
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;

  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel & Graphlet Count Kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-e editdistance : Set the Edit Distance computing : " << endl;
  cerr << "\t0 -- V1" << endl;
  cerr << "\t1 -- V2" << endl;
  cerr << "-N : Normalize Treelet distribution" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-o : Output gnuplot file. If not specified, gnuplot display" << endl;
  
}

bool always_true(double i)
{
  return true;
}

map<string, double,bool (*)(string, string)>* vit_to_weight(vector<string> vit[SIZE_SPECTRUM])
{
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  map<string, double,bool (*)(string, string)> * weights = 
    new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
  
  
  for(int k=0;k<SIZE_SPECTRUM;k++)
    {
      weights[k] = map<string, double, bool (*)(string, string)> (fn_pt); 
      for(unsigned int i=0;i<vit[k].size();i++)
	weights[k].insert(pair<string,double>(vit[k][i],1.0));
    }
  return weights;
}

vector<string> * getAllAsVIT (MoleculesDataset & dataset){
  dataset.computeLabeledSpectrums();
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  int nb_patterns = 0;
  //Initialisation de la vit list avec tous les treelets
  for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      treelet_spectrum::iterator it = distribution[treelet_type].begin();
      for(;it != distribution[treelet_type].end();it ++)
	{
	  vit_list[treelet_type].push_back(it->first);
	  nb_patterns ++;
	}
    }
  return vit_list;
}

void readOptions (int argc, char** argv, Options* options)
{
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  if (argc < 3)
    {
      usage(argv[0]);
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
  options->dataset_path = argv[1];
  options->dataset_file = argv[2];
  options->quietMode = false;
  options->output_file = NULL;
  options->normalize = false;
  options->spectrumKernel = InnerProductKernelType;

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
      else if (strncmp(argv[i], "-o", 2) == 0)
	{
	  options->output_file = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
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
	  usage(argv[0]);
	  exit(1);
	}
    }
}

int main (int argc, char** argv)
{	
  Options options;
  readOptions(argc, argv, &options);
  
  //Initializations
  ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
  GraphKernel* kgraph = NULL;
  GraphEditDistance* edit = NULL;
  if(options.editdistance == 0){
    edit = new MoleculeGraphEditDistance;
  }else{
    edit = new MoleculeGraphEditDistanceV2;
  }
  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.dataset_path, options.dataset_file);
  double * coeffs;
  //Kernel initialization
  switch(options.kernel)
    {
    case 0:
      cout << "Kernel : KMean" << endl;
      kgraph = new KMean(sbt, NULL);
      break;
			
    case 1:
      cout << "Kernel : KWMean" << endl;
      cout << "Sigma : " << options.sigma << endl;
      kgraph = new KWMean (sbt, NULL, options.sigma);
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
      kgraph = new LaplacianKernel (edit, dataset, options.sigma, options.regularization, options.lambda);
      break;
      
    case 3:
      cout << "Kernel : KEditDistance" << endl;
      kgraph = new KEditDistance (edit, options.sigma);
      break;	
    case 4:  
      cout << "Kernel : Unlabeled Treelet" << endl;
      dataset.computeSpectrums();
      dataset.computeGraphletCorrelation();
      if(!options.quietMode)
	dataset.printSpectrums();
      kgraph = new GraphletCountKernel(options.sigma,options.mu);
      coeffs = new double[SIZE_SPECTRUM];
      for(int i = 0;i < SIZE_SPECTRUM;i++)
	coeffs[i] = dataset.getCorrelationCoeff(i);
      ((GraphletCountKernel*)kgraph)->setCoeffs(coeffs,SIZE_SPECTRUM);
      break;
    case 5:
      cout << "Kernel : Random Walks Kernel" << endl;
      kgraph = NULL;//new RandomWalkKernel(new Kashima(),options.lambda,1.0);
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
      dataset.computeLabeledSpectrums();
      if(options.normalize)
	dataset.normalizeLabeledSpectrums();
      kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
      {
	vector<string> * vit_list = getAllAsVIT(dataset);
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
      }
	break;
    default:
      cerr << "Error : Invalid kernel" << endl;
      exit(1);
    }

  if(!options.quietMode){
    cout << "Dataset : " << options.dataset_path << "/" << options.dataset_file << endl;
    cout << "Number of molecules : " << dataset.size() << endl;
    cout << "Alpha : " << options.alpha << endl;
    cout << endl;
  }
  //End Initializations
  
  //Normalize for G'MIC script
  dataset.normalizeParams(0);
  dataset.computeGramMatrix(kgraph, false);
  KernelPCA * pca_projection = new KernelPCA(kgraph, dataset);
  double ** projection_result = pca_projection->projectTrainingSet(2);

  if(options.output_file == NULL)
    {
      //No output file given => Instant plotting
      double * x = new double[dataset.size()];
      double * y = new double[dataset.size()];
      // double * z = new double[dataset.size()];
      
      for(int i=0;i<dataset.size();i++){
	x[i] = projection_result[i][0];
	y[i] = projection_result[i][1];
	cout << i << " :" << x[i] << "," << y[i] << endl;
	// z[i] = projection_result[i][2];
      } 
      gnuplot_plot_once("Gram Matrix Projection", "points" ,"x","y",x,y,dataset.size()) ;
      delete [] x;
      delete [] y;
    }
  else
    {
      //DLM file generation
      fstream result_file;
      result_file.open(options.output_file,fstream::out | fstream::trunc);
      result_file.precision(DBL_DIG);
      for(int i = 0; i<dataset.size();i++)
	{
	  result_file << projection_result[i][0] << " " 
		      << projection_result[i][1] << " " 
		      << dataset.getParameter(i) << " " 
		      << dataset.getMoleculeFilename(i) << endl;
	}
      result_file.close();
    }

  delete projection_result;
  dataset.destroy();
  delete edit;
  delete kgraph;

  return 0;
}

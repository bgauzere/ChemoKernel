/*
 * @file moleculesClassification.cpp
 *
 * This program performs a SVM classification on molecules datasets.
 *
 * If no test set is given, then a leave-one-out cross-validation will be
 * performed on the training set.
 *
 * The library LibSVM (http://www.csie.ntu.edu.tw/~cjlin/libsvm/) is used for the
 * classification.
 *
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-10)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cfloat>
#include <sys/times.h>
#include <sys/time.h>

#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "TrailKernel.h"
#include "ShapeBagTrails.h"
#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "LaplacianKernelOriginal.h"
#include "GraphEditDistance.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"
#include "GraphletCountKernel.h"
#include "Kashima.h"
#include "KEditDistance.h"
#include "TreeletKernel.h"
#include "TanimotoKernel.h"
#include "MinMaxKernel.h"
#include <SVM.h>
#include "string_utils.h"
#include "KPrototypeEditDistance.h"
#include "CycleKernel.h"
#include "AugmentedCycleKernel.h"
#include "ContractedCycleKernel.h"
#include "TreeletCycleKernel.h"
#include "CombinedKernel.h"
#include "HorvathKernel.h"
#include "InterTreeletKernel.h"
#include "GaussianKernel.h"
#include "IntersectionKernel.h"
#include "InnerProductKernel.h"
#include "Kernel.h"
#include "CImg.h"
#include "TreeletEditDistance.h"
#include "MoleculeGraphEditDistanceMCS.h"
#include "TreeletType.h"
#include "TreeletIntersectionKernel.h"
#include "TreeletGaussianKernel.h"
#include "TreeletInnerProductKernel.h"
#include "TreeletBinaryKernel.h"
#include "TreeletRandomKernel.h"
#include "TreeletCompleteGaussianKernel.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

//XXX: Code a factoriser
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
vector<string> * getSizeTreeletsasVIT (MoleculesDataset & dataset, int taille){
  // int list_struct[14];
  // int nb_treelets = 0;
  // switch (taille){
  // case 1:
  //   list_struct = {0};
  //   nb_treelets = 1;    
  //   break;
  // case 2:
  //   list_struct = {0,1};
  //   nb_treelets = 2;    
  //   break;
  // case 3:
  //   list_struct = {0,1,2};
  //   nb_treelets = 3;    
  //   break;
  // case 4:
  //   list_struct = {0,1,2,3,6};
  //   nb_treelets = 5;    
  //   break;
  // case 5:
  //   list_struct = {0,1,2,3,4,6,7,8};
  //   nb_treelets = 8;    
  //   break;
  // case 7:
  //   list_struct = {0,1,2,3,4,5};
  //   nb_treelets = 6;    
  //   break;
  // default:
  //   list_struct = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
  //   nb_treelets = 14;    
  //   ;
  // }
  // dataset.computeLabeledSpectrums();
  // treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  // int nb_patterns = 0;
  // //Initialisation de la vit list avec tous les treelets
  // for(int i =0;i<nb_treelets;i++)
  //   {
  //     treelet_spectrum::iterator it = distribution[list_struct[i]].begin();
  //     for(;it != distribution[list_struct[i]].end();it ++)
  // 	{
  // 	  vit_list[list_struct[i]].push_back(it->first);
  // 	  nb_patterns ++;
  // 	}
  //   }
  return vit_list;
}


vector<string> * getAllAsVIT (MoleculesDataset & dataset){
  dataset.computeLabeledSpectrums();
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  int N = dataset.size();
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
vector<string> * getCyclesAllAsVIT (MoleculesDataset & dataset){
  treelet_spectrum * distribution = dataset.getCycleTreeletDistribution(&always_true);
  int N = dataset.size();
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

vector<string> * getAugmentedCyclesAllAsVIT (MoleculesDataset & dataset){
  treelet_spectrum * distribution = dataset.getAugmentedCycleTreeletDistribution(&always_true);
  int N = dataset.size();
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

vector<string> * getContractedCyclesAllAsVIT (MoleculesDataset & dataset){
  treelet_spectrum * distribution = dataset.getContractedCycleTreeletDistribution(&always_true);
  int N = dataset.size();
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

//Fin code a factoriser

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
  double c;
  bool useTestset;
  bool quietMode;
  int customSizeTrainset;
  bool normalize;
  KernelType spectrumKernel;
  KernelType spectrumKernelCycle;
  KernelType spectrumKernelContractedCycle;
  KPrototypeEditDistance::PrototypeSelectionType prototypeSelection;
  char * weights_file;
  char * grams_file;
  double cs,ci; //costs for substitution and insertion/deletion for edit distance
  double sigma_inter;
  char * distance_file;
  char * output_prefix;
  int size_max;
  int labelCycles;
};

void showUsage ()
{
  cerr << "Usage : moleculesClassification trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
	
  cerr << "-k kernel_type : Set the type of kernel function (default 0)" << endl;
  cerr << "\t0 -- Kernel Mean" << endl;
  cerr << "\t1 -- Kernel Weighted Mean" << endl;
  cerr << "\t2 -- Graph Laplacian Kernel" << endl;
  cerr << "\t3 -- Gaussian Kernel" << endl;
  cerr << "\t4 -- Unlabeled Treelet Kernel" << endl;
  cerr << "\t5 -- Random Walk Kernel" << endl;
  cerr << "\t6 -- Treelet Kernel" << endl;  
  cerr << "\t7 -- Tanimoto Kernel on Treelets" << endl;
  cerr << "\t8 -- MinMax Kernel on Treelets" << endl;
  cerr << "\t9 -- Test Graph Laplacian Kernel" << endl;
  cerr << "\t10 -- Edit Distance Prototype Selection" << endl;
  cerr << "\t11 -- Cycle Kernel" << endl;
  cerr << "\t12 -- Treelet + Cycle Kernel (-w to define weights, -l for treelet kernel weight and -s for cycle kernel)" << endl;
  cerr << "\t13 -- All Kernels (gaussian 1,2,3, inner, inter)" << endl;
  cerr << "\t14 -- Weighted Grams (use with -g)" << endl;
  cerr << "\t15 -- Horvath Kernel" << endl;

  cerr << "\t16 -- Inter Treelet Kernel - Approximate Edit Distance (No Regul)" << endl;
  cerr << "\t17 -- Inter Treelet Kernel - Laplacian Kernel" << endl;
  cerr << "\t18 -- Inter Treelet Kernel - Exact Connexe Edit Distance (Regul)" << endl;
  cerr << "\t21 -- Inter Treelet Kernel - avec distance_file" << endl;
  cerr << "\t19 -- Treelet Kernel on Contracted Cycle Hypergraph" << endl;

  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "-M kernel_type : Set the type of kernel between the two spectrums extracted from Contracted Cycle Hypergraph (CCH Treelet Kernel)" << endl;
  cerr << "-L kernel_type : Set the type of kernel between the two spectrums (Cycle Kernel Combo)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;
  cerr << "\t6 -- Inner Gaussian Kernel" << endl;


  cerr << "-P method : Set the method to select the prototypes" << endl;
  cerr << "-r nb proto : using -k 10, set the number of prototypes" << endl;

  cerr << "-t testset_path dataset_file : Set the test set to use (default none)" << endl;
  cerr << "-c C : Set the SVM parameter C (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-n : Size of the train dataset considered to predict the class of the molecule tested" << endl;
  cerr << "-N : Normalize Kernel" << endl;
  cerr << "-w : weights file" << endl;
  cerr << "-g : gram and weights file" << endl;
  cerr << "-ci : Node/Edge Insertion/Deletion cost" << endl;
  cerr << "-cs : Node/Edge substitution cost" << endl;
  cerr << "-sigma_inter : Sigma pour la similarité inter treelets" << endl;
  cerr << "-distance_file file : fichier encodant la distance entre treelets pour le train set" << endl;
  cerr << "-size_max  size : Taille maximale des treelets" << endl;
  cerr << "-lC Modifications of label"<< endl;
  cerr << "\t0 -- No modifications" << endl;
  cerr << "\t1 -- Number of relevant cycles add to the label of an atom" << endl;
  cerr << "\t2 -- Size of each relevant cycles add to the label of an atom" << endl;

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
  options->c = 1.0;
  options->useTestset = false;
  options->quietMode = false;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->customSizeTrainset = 0;
  options->normalize = false;
  options->spectrumKernel = InnerProductKernelType;  
  options->spectrumKernelContractedCycle = InnerProductKernelType;  
  options->spectrumKernelCycle = InnerProductKernelType;
  options->weights_file = NULL;
  options->distance_file = NULL;
  options->grams_file = NULL;
  options->output_prefix = NULL;
  options->ci = 7;
  options->cs = 1;
  options->sigma_inter = 1;
  options->size_max = -1;
  options->labelCycles=0;

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
      else if (strncmp(argv[i], "-lC", 3) == 0)
	{
	  options->labelCycles = atoi(argv[i+1]);
	  i+=2;
	} 
      else if (strcmp(argv[i], "-s") == 0)
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
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atof(argv[i+1]);
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
      else if (strncmp(argv[i], "-w", 2) == 0)
	{
	  options->weights_file = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-g", 2) == 0)
	{
	  options->grams_file = argv[i+1];
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
      else if (strncmp(argv[i], "-M", 2) == 0)
	{
	  options->spectrumKernelContractedCycle = (KernelType)(atoi(argv[i+1]));
	  i+=2;
	}
      else if (strncmp(argv[i], "-L", 2) == 0)
	{
	  options->spectrumKernelCycle = (KernelType)(atoi(argv[i+1]));
	  i+=2;
	}
      else if (strncmp(argv[i], "-P", 2) == 0)
	{
	  options->prototypeSelection = (KPrototypeEditDistance::PrototypeSelectionType)(atoi(argv[i+1]));
	  i+=2;
	}
      else if (strcmp(argv[i], "-ci") == 0)
	{
	  options->ci = atof(argv[i+1]);
	  i+=2;
	}
      else if (strcmp(argv[i], "-cs") == 0)
	{
	  options->cs = atof(argv[i+1]);
	  i+=2;
	}
      else if (strcmp(argv[i], "-sigma_inter") == 0)
	{
	  options->sigma_inter = atof(argv[i+1]);
	  i+=2;
	}
      else if (strcmp(argv[i], "-distance_file") == 0)
	{
	  options->distance_file = argv[i+1];
	  i+=2;
	}
      else if (strcmp(argv[i], "-output_prefix") == 0)
	{
	  options->output_prefix = argv[i+1];
      i+=2;
	}
      else if (strcmp(argv[i], "-size_max") == 0)
	{
	  options->size_max = atof(argv[i+1]);
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
	
  CycleKernel * kernel_c;
  ContractedCycleKernel * kernel_ccH;
  TreeletKernel * kernel_t;
  long diff;
  clock_t start,end;
  struct timeval s_start, s_end;

  long diff_total;
  clock_t start_total,end_total;
  struct timeval s_start_total, s_end_total;

  start_total = gettimeofday(&s_start_total,NULL);

  switch(options.kernel)
    {
    case 0:
      if(!options. quietMode)
	cout << "Kernel : KMean" << endl;
      kgraph = new KMean(sbt, NULL); //new Kashima()
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
      cout << "Kernel : KEditDistance" << endl;
      kgraph = new KEditDistance (edit, options.sigma);
      break;	
    case 4:
      cout << "Kernel : Unlabeled Treelet" << endl;
      trainset.computeSpectrums();
      trainset.printSpectrums();
      kgraph = new GraphletCountKernel();
      break;
    case 5:
      cout << "Kernel : Random Walks Kernel" << endl;
      kgraph = new RandomWalkKernel(new Kashima(),options.lambda,1.0);
      break;			
    case 6:  
      cout << "Kernel : Labeled Treelet" << endl;
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
	case CompleteGaussianKernelType:
	  cout << "CompleteGaussianKernelType";
	  break;
	case InnerGaussianKernelType:
	  cout << "InnerGaussianKernelType";
	  break;
	case PolynomialKernelType:
	  cout << "PolynomialKernelType";
	  break;
	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      trainset.computeLabeledSpectrums();
      if(options.normalize)
	trainset.normalizeLabeledSpectrums();
      if(options.useTestset){
	testset.computeLabeledSpectrums();
	if(options.normalize)
	  testset.normalizeLabeledSpectrums();
      }
      
      kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
      if(options.weights_file == NULL){
	((TreeletKernel*)kgraph)->selectTreelets(getAllAsVIT(trainset));
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(getAllAsVIT(trainset)));
      }else{
	vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
	bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	map<string, double,bool (*)(string, string)> * weights = 
	  new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	for(int k=0;k<SIZE_SPECTRUM;k++)
	  weights[k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	
	//lecture du fichier
	ifstream f;
	f.open(options.weights_file, ifstream::in );
	assert(f.good());
	string line;
	while(getline(f,line))
	  {
	    stringstream parse_line;
	    string code;
	    double weight;
	    
	    int pos_tiret = line.find("-");
	    int pos_sep_code_weight = line.find_last_of(" ");
	    
	    parse_line << line.substr(0,pos_tiret);
	    int treelet_type;
	    parse_line >> treelet_type;

	    code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	    // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	    // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	    weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	    weights[treelet_type].insert(pair<string, double>(code,weight));
	    vit_list[treelet_type].push_back(code);
	  }
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(weights);
      }
      if(options.size_max > 0)
	{
	  vector <string> * vit_list = getSizeTreeletsasVIT(trainset, options.size_max);
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      trainset.getTreeletDistribution(&always_true);
      cout << "Nb Treelets : " << trainset.getNbTreelets() << endl;
      testset.getTreeletDistribution(&always_true);
      cout << "Nb Treelets : " << testset.getNbTreelets() << endl;
      break;
    case 7:  
      cout << "Kernel : Tanimoto Kernel on Treelets" << endl;
      trainset.computeLabeledSpectrums();
      kgraph = new TanimotoKernel();
      break;
    case 8:  
      cout << "Kernel : MinMax Kernel on Treelets" << endl;
      trainset.computeLabeledSpectrums();
      kgraph = new MinMaxKernel();
      break;
    case 9:
      if(!options.quietMode){
	cout << "Kernel : Graph Laplacian Kernel (Test Mode)" << endl;
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
	kgraph = new LaplacianKernelOriginal (edit, trainset, testset, 
					      options.sigma, options.regularization, options.lambda);
      else
	kgraph = new LaplacianKernelOriginal (edit, trainset, options.sigma, 
					      options.regularization, options.lambda);
      break;
    case 10:
      cout << "Kernel : Prototype Edit Distance" << endl;
      if (options.useTestset)
	kgraph = new KPrototypeEditDistance(trainset, testset, options.regularization, options.prototypeSelection);
      
      else
	kgraph = new KPrototypeEditDistance(trainset, options.regularization,options.prototypeSelection);
      break;
    case 11:
      cout << "Kernel : Cycle Treelet" << endl;
      cout << "Spectrum Kernel : ";
      switch(options.spectrumKernelCycle)
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
	case CompleteGaussianKernelType:
	  cout << "CompleteGaussianKernelType";
	  break;
	case InnerGaussianKernelType:
	  cout << "InnerGaussianKernelType";
	  break;
	case PolynomialKernelType:
	  cout << "PolynomialKernelType";
	  break;
	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      trainset.computeRelevantCycles();
      start = gettimeofday(&s_start,NULL);
      trainset.computeCyclesSpectrums();
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Cycle Spectrum Computation Time : " << diff << endl;
      
      if(options.useTestset){
	testset.computeRelevantCycles();
	testset.computeCyclesSpectrums();
      }
      kgraph = new CycleKernel(options.spectrumKernelCycle,options.sigma);
      
      if(options.weights_file == NULL){
	((CycleKernel*)kgraph)->selectTreelets(getCyclesAllAsVIT(trainset));
	((CycleKernel*)kgraph)->weightTreelets(vit_to_weight(getCyclesAllAsVIT(trainset)));
      }else{
	vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
	bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	map<string, double,bool (*)(string, string)> * weights = 
	  new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	for(int k=0;k<SIZE_SPECTRUM;k++)
	  weights[k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	
	//lecture du fichier
	ifstream f;
	f.open(options.weights_file, ifstream::in );
	assert(f.good());
	string line;
	while(getline(f,line))
	  {
	    stringstream parse_line;
	    string code;
	    double weight;
	    
	    int pos_tiret = line.find("-");
	    int pos_sep_code_weight = line.find_last_of(" ");
	    
	    parse_line << line.substr(0,pos_tiret);
	    int treelet_type;
	    parse_line >> treelet_type;

	    code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	    // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	    // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	    weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	    weights[treelet_type].insert(pair<string, double>(code,weight));
	    vit_list[treelet_type].push_back(code);
	  }
	((CycleKernel*)kgraph)->selectTreelets(vit_list);
	((CycleKernel*)kgraph)->weightTreelets(weights);
      }
      trainset.getCycleTreeletDistribution(&always_true);
      cout << "Nb Cycle Treelets : " << trainset.getNbCyclesTreelets() << endl;
      testset.getCycleTreeletDistribution(&always_true);
      cout << "Nb Cycle Treelets : " << testset.getNbCyclesTreelets() << endl;
      break;
    case 12:
      trainset.computeLabeledSpectrums();
      trainset.computeRelevantCycles();
      trainset.computeCyclesSpectrums();
      if(options.useTestset){
	testset.computeLabeledSpectrums();
	testset.computeRelevantCycles();
	testset.computeCyclesSpectrums();
      }
      kernel_t = new TreeletKernel(options.spectrumKernel,options.sigma);
      kernel_c = new CycleKernel(options.spectrumKernelCycle,options.sigma);
      if(options.weights_file == NULL){
	kernel_c->selectTreelets(getCyclesAllAsVIT(trainset));
	kernel_c->weightTreelets(vit_to_weight(getCyclesAllAsVIT(trainset)));
	kernel_t->selectTreelets(getAllAsVIT(trainset));
	kernel_t->weightTreelets(vit_to_weight(getAllAsVIT(trainset)));
      }else{
	//XXX: A factoriser
	vector<vector<string> *> vit_lists(2);
	bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	vector<map<string, double,bool (*)(string, string)> *> weights(2);
	for(int i=0;i<2;i++){
	  vit_lists[i] = new vector<string>[SIZE_SPECTRUM];
	  weights[i] = new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	  for(int k=0;k<SIZE_SPECTRUM;k++)
	    weights[i][k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	}
	//lecture du fichier
	int n_treelet =0;
	int n_kernel = 0;
	ifstream f;
	f.open(options.weights_file, ifstream::in );
	assert(f.good());
	string line;
	while(getline(f,line))
	  {
	    n_treelet ++;
	    if(line[0] == '-'){
	      n_treelet = 0;
	      n_kernel++;
	      cerr << "Treelet -> Cycle" << endl;
	    }else{
	      stringstream parse_line;
	      string code;
	      double weight;
	      
	      int pos_tiret = line.find("-");
	      int pos_sep_code_weight = line.find_last_of(" ");
	    
	      parse_line << line.substr(0,pos_tiret);
	      int treelet_type;
	      parse_line >> treelet_type;

	      code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	      // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	      // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	      weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	      weights[n_kernel][treelet_type].insert(pair<string, double>(code,weight));
	      vit_lists[n_kernel][treelet_type].push_back(code);
	    }
	  }
	kernel_t->selectTreelets(vit_lists[0]);
	kernel_t->weightTreelets(weights[0]);
	kernel_c->selectTreelets(vit_lists[1]);
	kernel_c->weightTreelets(weights[1]);
      }
      // trainset.computeGramMatrix(kernel_t, false);
      // trainset.showGramMatrixMatlab("treelets.mat");
      // trainset.computeGramMatrix(kernel_c, false);
      // trainset.showGramMatrixMatlab("cycles.mat");
      // exit(0);
      kgraph = new TreeletCycleKernel(kernel_t,kernel_c, options.lambda, 1-options.lambda);
      break;
    case 13:
      cout << "All Kernels "<< endl;
      trainset.computeLabeledSpectrums();
      if(options.normalize)
	trainset.normalizeLabeledSpectrums();
      if(options.useTestset){
	testset.computeLabeledSpectrums();
	if(options.normalize)
	  testset.normalizeLabeledSpectrums();
      }
      {
	//Construction des 5 kernels
	vector <GraphKernel *> kernels(5);
	kernels[0] = new TreeletKernel(IntersectionKernelType);
	kernels[1] = new TreeletKernel(GaussianKernelType,0.5);
	kernels[2] = new TreeletKernel(GaussianKernelType,1);
	kernels[3] = new TreeletKernel(GaussianKernelType,1.5);
	kernels[4] = new TreeletKernel(InnerProductKernelType);

	vector<double> weights(5,1.0);
	if(options.weights_file == NULL){
	  for(int i=0;i<kernels.size();i++){
	    ((TreeletKernel*)kernels[i])->selectTreelets(getAllAsVIT(trainset));
	    ((TreeletKernel*)kernels[i])->weightTreelets(vit_to_weight(getAllAsVIT(trainset)));
	  }
	}else{
	  vector<vector<string> *> vit_lists(kernels.size());
	  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	  vector<map<string, double,bool (*)(string, string)> *> weights(kernels.size());
	  for(int i=0;i<kernels.size();i++){
	    vit_lists[i] = new vector<string>[SIZE_SPECTRUM];
	    weights[i] = new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	    for(int k=0;k<SIZE_SPECTRUM;k++)
	      weights[i][k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	  }
	
	  //lecture du fichier
	  int n_treelet =0;
	  int n_kernel = 0;
	  ifstream f;
	  f.open(options.weights_file, ifstream::in );
	  assert(f.good());
	  string line;
	  while(getline(f,line))
	    {
	      n_treelet ++;
	      if(n_treelet == 4676){
		n_treelet =0;
		n_kernel++;
	      }
	      stringstream parse_line;
	      string code;
	      double weight;
	    
	      int pos_tiret = line.find("-");
	      int pos_sep_code_weight = line.find_last_of(" ");
	    
	      parse_line << line.substr(0,pos_tiret);
	      int treelet_type;
	      parse_line >> treelet_type;

	      code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	      // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	      // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	      weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	      weights[n_kernel][treelet_type].insert(pair<string, double>(code,weight));
	      vit_lists[n_kernel][treelet_type].push_back(code);
	    }
	  for(unsigned int i=0;i<kernels.size();i++){
	    ((TreeletKernel*)kernels[i])->selectTreelets(vit_lists[i]);
	    ((TreeletKernel*)kernels[i])->weightTreelets(weights[i]);
	  }
	}
	kgraph = new CombinedKernel(kernels,weights);
	break;
      }    
    case 14:
      cout << "Not implemented yet."<< endl;
      exit(0);
      break;
    case 15:
      cout << "Horvath Kernel"<< endl;
      trainset.computeRelevantCycles();
      
      start = gettimeofday(&s_start,NULL);
      trainset.computeSimpleCycles(options.sigma);
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Simple Cycles computation Time : " << diff << endl;
      if(options.useTestset){
	testset.computeRelevantCycles();
	testset.computeSimpleCycles(options.sigma);
      }
      kgraph = new HorvathKernel();
      break;
    case 16:
      {
	trainset.computeLabeledSpectrums();
	if(options.useTestset){
	  testset.computeLabeledSpectrums();
	}
	treelet_spectrum * distrib = trainset.getTreeletDistribution(&always_true);
	int nb_treelets =  trainset.getNbTreelets();
	cout << nb_treelets<< endl;
	trainset.computeSpecialVectors(&distrib, nb_treelets);
	if(options.useTestset){
	  testset.computeSpecialVectors(&distrib, nb_treelets);
	}
	cout << "Special Vectors computed " << endl;
	
	Kernel * kinter_graphs = NULL;
	switch(options.spectrumKernel){
	case IntersectionKernelType:
	  kinter_graphs = new IntersectionKernel();
	  break;
	case InnerProductKernelType:
	  kinter_graphs = new InnerProductKernel();
	  break;
	case GaussianKernelType:
	  kinter_graphs = new GaussianKernel(options.sigma);
	  break;
	default:
	  exit(EXIT_FAILURE);
	}
	cout << "Edit Distance Initialization ..." << endl;
	clock_t itk_start,itk_end;
	struct timeval itk_s_start, itk_s_end;
	itk_start = gettimeofday(&itk_s_start,NULL);
	GraphEditDistance * ed_treelets = new MoleculeGraphEditDistanceMCS(options.cs, options.ci);
	cout << "Edit Distance Initialized" << endl;
	GraphKernel * kinter_treelets = new KEditDistance (ed_treelets, options.sigma_inter);
	itk_end = gettimeofday(&itk_s_end,NULL);//times(&s_end);
	long itk_time = (1000000*itk_s_end.tv_sec + itk_s_end.tv_usec) - 
	(1000000*itk_s_start.tv_sec + itk_s_start.tv_usec);	
	cout << "Inter Treelet Kernel initialized in " << itk_time << "micro seconds."<< endl;
	cout << "Edit distance Computed"<< endl;
	kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
					&distrib,
					nb_treelets);
      }
      break;
    case 17:      
      {
	trainset.computeLabeledSpectrums();
	if(options.useTestset){
	  testset.computeLabeledSpectrums();
	}
	treelet_spectrum * distrib = trainset.getTreeletDistribution(&always_true);
	int nb_treelets =  trainset.getNbTreelets();
	cout << nb_treelets<< endl;
	trainset.computeSpecialVectors(&distrib, nb_treelets);
	if(options.useTestset){
	  testset.computeSpecialVectors(&distrib, nb_treelets);
	}
	cout << "Special Vectors computed " << endl;
	clock_t itk_start,itk_end;
	struct timeval itk_s_start, itk_s_end;
	itk_start = gettimeofday(&itk_s_start,NULL);

	Dataset dataset_treelets;
	Collection ** col_treelets = new Collection* [nb_treelets];
	int cur_treelet = 0;
	for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
	  {
	    treelet_spectrum::iterator it = distrib[i].begin();
	    for(;it != distrib[i].end();it++,cur_treelet++){
	      //cout << i << "," << it->first << endl;
	      Collection * tmp = TreeletEnumerator::TreeletToCollection(i,it->first);
	      dataset_treelets.add(tmp,1.0,NULL);
	      col_treelets[cur_treelet] = tmp;
	    }
	  }
	cout << "Collections added" << endl;
	GraphEditDistance * ed_treelets = new MoleculeGraphEditDistanceMCS(options.cs, options.ci);
	GraphKernel * kinter_treelets = new LaplacianKernelOriginal (ed_treelets, dataset_treelets,
								     options.sigma_inter,
								     options.regularization,
								     options.lambda);
	cout << "Laplacian Computed" << endl;
	Kernel * kinter_graphs = NULL;
	switch(options.spectrumKernel){
	case IntersectionKernelType:
	  kinter_graphs = new IntersectionKernel();
	  break;
	case InnerProductKernelType:
	  kinter_graphs = new InnerProductKernel();
	  break;
	case GaussianKernelType:
	  kinter_graphs = new GaussianKernel(options.sigma);
	  break;
	default:
	  exit(EXIT_FAILURE);
	}
	cout << "TreeletKernel Initialized" << endl;
	kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
					col_treelets,
					nb_treelets);
	itk_end = gettimeofday(&itk_s_end,NULL);//times(&s_end);
	long itk_time = (1000000*itk_s_end.tv_sec + itk_s_end.tv_usec) - 
	(1000000*itk_s_start.tv_sec + itk_s_start.tv_usec);	
	cout << "Inter Treelet Kernel initialized in " << itk_time << "micro seconds."<< endl;
      }
      break;
    case 18:      
      {
	trainset.computeLabeledSpectrums();
	if(options.useTestset){
	  testset.computeLabeledSpectrums();
	}
	treelet_spectrum * distrib = trainset.getTreeletDistribution(&always_true);
	int nb_treelets =  trainset.getNbTreelets();
	cout << nb_treelets<< endl;
	trainset.computeSpecialVectors(&distrib, nb_treelets);
	if(options.useTestset){
	  testset.computeSpecialVectors(&distrib, nb_treelets);
	}
	cout << "Special Vectors computed " << endl;
	Kernel * kinter_graphs = NULL;
	switch(options.spectrumKernel){
	case IntersectionKernelType:
	  kinter_graphs = new IntersectionKernel();
	  break;
	case InnerProductKernelType:
	  kinter_graphs = new InnerProductKernel();
	  break;
	case GaussianKernelType:
	  kinter_graphs = new GaussianKernel(options.sigma);
	  break;
	default:
	  exit(EXIT_FAILURE);
	}
	cout << "Kernel Initialization ..." << endl;

	kgraph = new InterTreeletKernel(NULL,kinter_graphs,
					&distrib,
					nb_treelets);
	cout << "Kernel Initialized ..." << endl;
	cout << "Exact Edit Distance Computation ..." << endl;
	int * type_treelet = new int[nb_treelets];
	string* code_treelet = new string[nb_treelets];
	int cur_treelet = 0;
	for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
	  {
	    treelet_spectrum::iterator it = distrib[i].begin();
	    for(;it != distrib[i].end();it++){
	      //cout << i << "," << it->first << endl;
	      type_treelet[cur_treelet] = i;
	      code_treelet[cur_treelet++] = it->first;
	      // cout << i << " :" << MoleculeGraph::translateTreeletCode(it->first) << endl;
	    }
	  }
        CImg<double> sim_treelets(nb_treelets,nb_treelets,1,1,0.0);
	// for (unsigned int i=0; i<nb_treelets; ++i)
	//   for(int j=0;j<nb_treelets;j++)
	//     sim_treelets(i,j) = 1;
	cout << "Treelet Edit Distance Initialization ..." << endl;
	clock_t itk_start,itk_end;
	struct timeval itk_s_start, itk_s_end;
	itk_start = gettimeofday(&itk_s_start,NULL);
	TreeletEditDistance ed(options.cs,options.ci);
	TreeletEditDistance::Init();
	cout << "Treelet Edit Distance Initialized ..." << endl;
	cout << "Treelet Edit Distance Computation ..." << endl;
	for(int i=0;i<nb_treelets;i++)
	  for(int j=i;j<nb_treelets;j++){
	    sim_treelets(j,i) = sim_treelets(i,j) = exp(- pow(ed(type_treelet[i],code_treelet[i],
								 type_treelet[j],code_treelet[j]),2)/(options.sigma_inter));
	    if(sim_treelets(i,j) < DBL_EPSILON)
	      sim_treelets(i,j) = 0.0;
	  }
	
	((InterTreeletKernel*)kgraph)->setSimTreelets(sim_treelets);
	itk_end = gettimeofday(&itk_s_end,NULL);//times(&s_end);
	long itk_time = (1000000*itk_s_end.tv_sec + itk_s_end.tv_usec) - 
	  (1000000*itk_s_start.tv_sec + itk_s_start.tv_usec);	
	cout << "Inter Treelet Kernel initialized in " << itk_time << " micro seconds."<< endl;
	cout << "Treelet Edit Distance Computed ..." << endl;
	CImg<double> eigvals;
	CImg<double> eigvects;
	cout << "Calcul eigen values ..." << endl;
	sim_treelets.symmetric_eigen(eigvals,eigvects);
	cout << "regul ..." << endl;
	// for(int i=0;i<nb_treelets;i++)
	//   for(int j=0;j<nb_treelets;j++)
	//     assert(sim_treelets(i,j) == sim_treelets(j,i));
	if(eigvals[nb_treelets-1] < 0)
	  for (unsigned int i=0; i<nb_treelets; ++i)
	    sim_treelets(i,i) -= eigvals[nb_treelets-1];
	//sim_treelets.symmetric_eigen(eigvals,eigvects);
	//sim_treelets.display();
      }
      break;
    case 19:  
      cout << "Kernel : Labeled Treelet on Contracted Cycle Hypergraph" << endl;
      cout << "Spectrum Kernel : ";
      switch(options.spectrumKernelContractedCycle)
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
	case CompleteGaussianKernelType:
	  cout << "CompleteGaussianKernelType";
	  break;
	case InnerGaussianKernelType:
	  cout << "InnerGaussianKernelType";
	  break;
	case PolynomialKernelType:
	  cout << "PolynomialKernelType";
	  break;
	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      start = gettimeofday(&s_start,NULL);
      trainset.computeCCHypergraph();
      trainset.computeContractedCycleSpectrums();
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Contracted Cycle Computation Time : " << diff << endl;
      
      if(options.useTestset){
	testset.computeCCHypergraph();
	testset.computeContractedCycleSpectrums();
      }
      kgraph = new ContractedCycleKernel(options.spectrumKernelContractedCycle,options.sigma);
      if(options.weights_file == NULL){
	((ContractedCycleKernel*)kgraph)->selectTreelets(getContractedCyclesAllAsVIT(trainset));
	((ContractedCycleKernel*)kgraph)->weightTreelets(vit_to_weight(getContractedCyclesAllAsVIT(trainset)));
      }else{
	vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
	bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	map<string, double,bool (*)(string, string)> * weights = 
	  new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	for(int k=0;k<SIZE_SPECTRUM;k++)
	  weights[k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	
	//lecture du fichier
	ifstream f;
	f.open(options.weights_file, ifstream::in );
	assert(f.good());
	string line;
	while(getline(f,line))
	  {
	    stringstream parse_line;
	    string code;
	    double weight;
	    
	    int pos_tiret = line.find("-");
	    int pos_sep_code_weight = line.find_last_of(" ");
	    
	    parse_line << line.substr(0,pos_tiret);
	    int treelet_type;
	    parse_line >> treelet_type;

	    code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	    // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	    // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	    weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	    weights[treelet_type].insert(pair<string, double>(code,weight));
	    vit_list[treelet_type].push_back(code);
	  }
	((ContractedCycleKernel*)kgraph)->selectTreelets(vit_list);
	((ContractedCycleKernel*)kgraph)->weightTreelets(weights);

      }
      break; 
    case 20:
      trainset.computeLabeledSpectrums();
      trainset.computeCCHypergraph();
      trainset.computeContractedCycleSpectrums();
      if(options.useTestset){
	testset.computeLabeledSpectrums();
	testset.computeCCHypergraph();
	testset.computeContractedCycleSpectrums();
      }
      kernel_t = new TreeletKernel(options.spectrumKernel,options.sigma);
      kernel_ccH = new ContractedCycleKernel(options.spectrumKernelContractedCycle,options.sigma);
      if(options.weights_file == NULL){
	kernel_ccH->selectTreelets(getContractedCyclesAllAsVIT(trainset));
	kernel_ccH->weightTreelets(vit_to_weight(getContractedCyclesAllAsVIT(trainset)));
	kernel_t->selectTreelets(getAllAsVIT(trainset));
	kernel_t->weightTreelets(vit_to_weight(getAllAsVIT(trainset)));
      }else{
	//XXX: A factoriser
	vector<vector<string> *> vit_lists(2);
	bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	vector<map<string, double,bool (*)(string, string)> *> weights(2);
	for(int i=0;i<2;i++){
	  vit_lists[i] = new vector<string>[SIZE_SPECTRUM];
	  weights[i] = new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	  for(int k=0;k<SIZE_SPECTRUM;k++)
	    weights[i][k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	}
	//lecture du fichier
	int n_treelet =0;
	int n_kernel = 0;
	ifstream f;
	f.open(options.weights_file, ifstream::in );
	assert(f.good());
	string line;
	while(getline(f,line))
	  {
	    n_treelet ++;
	    if(line[0] == '-'){
	      n_treelet = 0;
	      n_kernel++;
	      cerr << "Treelet -> Contracted Cycle Hypergraph" << endl;
	    }else{
	      stringstream parse_line;
	      string code;
	      double weight;
	      
	      int pos_tiret = line.find("-");
	      int pos_sep_code_weight = line.find_last_of(" ");
	    
	      parse_line << line.substr(0,pos_tiret);
	      int treelet_type;
	      parse_line >> treelet_type;

	      code = line.substr(pos_tiret+1,pos_sep_code_weight-pos_tiret-1);
	      // parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	      // cout << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1) << endl;
	      weight = atof(line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1).c_str());
	      weights[n_kernel][treelet_type].insert(pair<string, double>(code,weight));
	      vit_lists[n_kernel][treelet_type].push_back(code);
	    }
	  }
	kernel_t->selectTreelets(vit_lists[0]);
	kernel_t->weightTreelets(weights[0]);
	kernel_ccH->selectTreelets(vit_lists[1]);
	kernel_ccH->weightTreelets(weights[1]);
      }
      // trainset.computeGramMatrix(kernel_t, false);
      // trainset.showGramMatrixMatlab("treelets.mat");
      // trainset.computeGramMatrix(kernel_ccH, false);
      // trainset.showGramMatrixMatlab("cycles.mat");
      // exit(0);
      
      {
	vector<GraphKernel *> kernels;
	vector<double> weights;
	kernels.push_back(kernel_t);
	kernels.push_back(kernel_ccH);
	weights.push_back(options.lambda);
	weights.push_back(1-options.lambda);
	kgraph = new CombinedKernel(kernels,weights);
      }
      break;
    case 21:
      {
	trainset.computeLabeledSpectrums();
	if(options.useTestset){
	  testset.computeLabeledSpectrums();
	}
	treelet_spectrum * distrib = trainset.getTreeletDistribution(&always_true);
	int nb_treelets =  trainset.getNbTreelets();
	cout << nb_treelets<< endl;
	trainset.computeSpecialVectors(&distrib, nb_treelets);
	if(options.useTestset){
	  testset.computeSpecialVectors(&distrib, nb_treelets);
	}
	cout << "Special Vectors computed " << endl;
	Kernel * kinter_graphs = NULL;
	switch(options.spectrumKernel){
	case IntersectionKernelType:
	  kinter_graphs = new IntersectionKernel();
	  break;
	case InnerProductKernelType:
	  kinter_graphs = new InnerProductKernel();
	  break;
	case GaussianKernelType:
	  kinter_graphs = new GaussianKernel(options.sigma);
	  break;
	default:
	  exit(EXIT_FAILURE);
	}
	cout << "Kernel Initialization ..." << endl;

	kgraph = new InterTreeletKernel(NULL,kinter_graphs,
					&distrib,
					nb_treelets);
	cout << "Kernel Initialized ..." << endl;
	CImg<double> sim_treelets(nb_treelets,nb_treelets,1,1,0.0);
	ifstream file(options.distance_file,ios::in);
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++){
	    int d;
	    file >> d;
	    sim_treelets(i,j) = exp(- (d*d)/(options.sigma_inter));
	    if(sim_treelets(i,j) < DBL_EPSILON)
	      sim_treelets(i,j) = 0.0;
	  }
	((InterTreeletKernel*)kgraph)->setSimTreelets(sim_treelets);
      break;
      }
    case 22:
      {
      cout << "Kernel : Labeled Treelet Augmented Cycles" << endl;
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
	case CompleteGaussianKernelType:
	  cout << "CompleteGaussianKernelType";
	  break;
	case InnerGaussianKernelType:
	  cout << "InnerGaussianKernelType";
	  break;
	case PolynomialKernelType:
	  cout << "PolynomialKernelType";
	  break;
	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      trainset.computeCCHypergraph();
      trainset.computeAugmentedCycles();
      trainset.computeAugmentedCycleSpectrums();
      if(options.useTestset){
	testset.computeCCHypergraph();
	testset.computeAugmentedCycles();
	testset.computeAugmentedCycleSpectrums();
      }
      
      kgraph = new AugmentedCycleKernel(options.spectrumKernel,options.sigma);
      ((AugmentedCycleKernel*)kgraph)->selectTreelets(getAugmentedCyclesAllAsVIT(trainset));
      ((AugmentedCycleKernel*)kgraph)->weightTreelets(vit_to_weight(getAugmentedCyclesAllAsVIT(trainset)));
      break;
    }
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
  long diff_gram;
  start = gettimeofday(&s_start,NULL);
  trainset.computeGramMatrix(kgraph, false);
  end_total = gettimeofday(&s_end_total,NULL);//times(&s_end);
  diff_total = (1000000*s_end_total.tv_sec + s_end_total.tv_usec) - 
    (1000000*s_start_total.tv_sec + s_start_total.tv_usec);
  cout.precision(DBL_DIG);
  cout << "Learning Time : " << diff_total << endl;
  cout << "Gram Matrix Computation Time : " << diff_gram << endl;
  cout << "Total Time : " << diff_gram + diff << endl;

  end = gettimeofday(&s_end,NULL);//times(&s_end);
  diff_gram = (1000000*s_end.tv_sec + s_end.tv_usec) - 
    (1000000*s_start.tv_sec + s_start.tv_usec);
  cout.precision(DBL_DIG);
  cout << "Gram Matrix Computation Time : " << diff_gram << endl;
  cout << "Total Time : " << diff_gram + diff << endl;
  trainset.showGramMatrixMatlab("dataset.mat");
  

  trainset.showGramMatrixRaw();
  //if(!options.quietMode)
  // cout << "Gram Matrix Positive definite ? ! " << trainset.isGramMatrixPD() << endl;
  // cout << "Gram Matrix Symmetric ? ! " << trainset.isSymmetric() << endl;
  double accuracy = 0.0;
  double fp,tp,fn,tn;
  fp = tp = fn = tn = 0;
  unsigned int N;
  double accuracyT = 0.0;
  double fpT,tpT,fnT,tnT;
  fpT = tpT = fnT = tnT = 0;
  unsigned int NT=trainset.size();
	
  if (options.useTestset)
    {
      testset.computeLabeledSpectrums();
      SVM svm (trainset, kgraph, SVM::C_SVC, options.c);
      N = testset.size();
      for (unsigned int i=0; i<N; ++i)
	{
	  double c = testset.getParameter(i);
	  double x = svm.predict(testset.getCollection(i));
	  if(!options.quietMode)
	    cout << "Molecule " << i << " : " << c << " --> " << x << endl;
	  if (c == x)
	    ++accuracy;
	
	  if ((c == 1) && (x == 1))
	    tp ++;
	  if ((c == 1) && (x == 0))
	    fn ++;
	  if ((c == 0) && (x == 1))
	    fp ++;
	  if ((c == 0) && (x == 0))
	    tn ++;	
	}
       for (unsigned int i=0; i<NT; ++i)
	{
	  double c = trainset.getParameter(i);
	  double x = svm.predict(trainset.getCollection(i));
	  if(!options.quietMode)
	    cout << "Train Molecule " << i << " : " << c << " --> " << x << endl;
	  if (c == x)
	    ++accuracyT;
	
	  if ((c == 1) && (x == 1))
	    tpT ++;
	  if ((c == 1) && (x == 0))
	    fnT ++;
	  if ((c == 0) && (x == 1))
	    fpT ++;
	  if ((c == 0) && (x == 0))
	    tnT ++;	
	}      
      // Destruction of the test set
      testset.destroy();
    }
  else // Leave-one-out cross-validation
    {
      N = trainset.size();
      int nbToSuppress = (N - options.customSizeTrainset) + 1;

      double * c = new double[nbToSuppress];
      Collection** col = new Collection*[nbToSuppress];
      
      for (unsigned int i=0; i<N; ++i)
	{
	  //On supprime le premier
	  c[nbToSuppress - 1] = trainset.getParameter(0);
	  col[nbToSuppress - 1] = trainset.getCollection(0);
	  trainset.delete_first();
	  
	  //Reajustement taille base apprentissage
	  for(int k = (nbToSuppress - 2); k >= 0 ; --k)
	    {
	      c[k] = trainset.getParameter(trainset.size() -1);
	      col[k] = trainset.getCollection(trainset.size() -1);
	      trainset.delete_last();
	    }

	  SVM svm (trainset, kgraph, SVM::C_SVC, options.c);
	  
	  double x = svm.predict(col[0]);
	
	  if(!options.quietMode)
	    cout << "Molecule " << i << " : " << c[0] << " --> " << x;
			
	  if (c[0] == x)
	    ++accuracy;
	  else
	    if(!options.quietMode)
	      cout << "****" ;
	  if(!options.quietMode)
	    cout << endl;
	  
	  if ((c[0] == 1) && (x == 1))
	    tp ++;
	  if ((c[0] == 1) && (x == 0))
	    fn ++;
	  if ((c[0] == 0) && (x == 1))
	    fp ++;
	  if ((c[0] == 0) && (x == 0))
	    tn ++;
			   
	  //retablissement taille base apprentissage
	  for(int k = 0; k < nbToSuppress; ++k)
	    {
	      //On le rajoute a la fin
	      trainset.add (col[k], c[k], kgraph);
	    }
	}
    }

  if(!options.quietMode)
    cout << endl;

  if(options.quietMode)
    cout << options.customSizeTrainset << " " << 100*accuracy/N << endl;
  else
    {
      cout << "Accuracy : " << 100*accuracy/N << "% (" << accuracy << " / " << N << ")" << endl;
      double precision = tp/(tp+fp);
      double recall = tp/(tp+fn);
      double f_score = 2*(precision*recall) / (precision+recall);
      cout << precision << "," << recall << "," << fp << "," << tn << ","<< tp << "," << fn << endl;
      cout << "F Score : " << f_score << endl;

      cout << "Accuracy Train : " << 100*accuracyT/NT << "% (" << accuracyT << " / " << NT << ")" << endl;
      double precisionT = tpT/(tpT+fpT);
      double recallT = tpT/(tpT+fnT);
      double f_scoreT = 2*(precisionT*recallT) / (precisionT+recallT);
      cout << precisionT << "," << recallT << "," << fpT << "," << tnT << ","<< tpT << "," << fnT << endl;
      cout << "F Score : " << f_scoreT << endl;

    }
  
  // Destruction of the training set
	
  trainset.destroy();
	
  delete edit;
  delete kgraph;
	
  return 0;
}


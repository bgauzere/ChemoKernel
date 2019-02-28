#include <iostream>
#include <algorithm>
#include <cfloat>
#include <sstream>
#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include <pandore.h>
#include "TreeletKernel.h"
#include "MoleculesDataset.h"
#include <libgen.h>
#include "utils.h"
#include "WeightedRegression.h"
#include "CImg.h"
#include "string_utils.h"
#include "KernelRidge.h"
#include <fstream>
#include "TreeletType.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;
#define LIMIT_LOW_FREQUENCY 140
vector<string> * getAcyclicFromFile(string file)
{
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  ifstream vit_file (file.c_str());
  string line;
  while (getline(vit_file,line))
    {
      stringstream parse;
      parse << line.substr(0,line.find("-"));
      int treelet_type;
      parse >> treelet_type;

      string code_line = line.substr(line.find("-")+1, line.length() - line.find("-"));
      vit_list[treelet_type].push_back(code_line);
    }
    return vit_list;
}

vector<string> * getAcyclicVIT_2()
{
  /*Liste des VIT déterminés par le complete kernel sur le dataset sans les similaires*/
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];  
  
  vit_list[0].push_back("C");
  vit_list[0].push_back("S");

  vit_list[1].push_back("C1O");
  vit_list[1].push_back("C1C");
  vit_list[1].push_back("S1S");

  vit_list[2].push_back("O1C1O");
  vit_list[2].push_back("C1O1O");
  vit_list[2].push_back("S1C1S");

  vit_list[3].push_back("C1C1C1C");
  vit_list[3].push_back("C1O1O1C");

  vit_list[4].push_back("C1S1C1C1S");
  vit_list[4].push_back("O1C1C1C1O");
  vit_list[4].push_back("C1C1C1C1S");
  vit_list[4].push_back("C1S1C1S1C");
  

  vit_list[5].push_back("C1O1C1C1C1O");
  vit_list[5].push_back("C1C1C1C1C1C");
  vit_list[5].push_back("C1C1C1O1C1O");
  vit_list[5].push_back("C1C1O1O1C1C");
  vit_list[5].push_back("C1S1C1C1S1C");
  vit_list[5].push_back("C1C1C1O1C1C");
  vit_list[5].push_back("C1C1C1C1C1S");
  vit_list[5].push_back("O1C1C1C1C1O");
  vit_list[5].push_back("C1C1C1C1S1S");
  vit_list[5].push_back("C1C1C1S1C1S");
  
  vit_list[6].push_back("C1C1C1S");
  vit_list[6].push_back("C1C1C1O");
  vit_list[6].push_back("C1C1C1C");

  vit_list[7].push_back(" S1C1C1C1C");

  vit_list[8].push_back("C1C1C1O1O");

  vit_list[9].push_back("C1C1C1C1S1S");  

  vit_list[10].push_back("C1O O1C1C1C1C");

  vit_list[11].push_back("C1O1C1C1C1O");  
  vit_list[11].push_back("O1O1C1C1C1C");

  vit_list[12].push_back("1 C1C1C C1C1O");

  return vit_list;
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

bool comp (pair<string, double>  l, pair<string, double>  r)
{ 
  return (fabs(l.second) > fabs(r.second)); 
}

bool always_true(double i)
{
  i = i;//compiler warning tip
  return true;
}

struct Options
{
  char* trainset_file;
  char* testset_file;
  
  double sigma;
  double lambda;
  double mu;
  
  bool quietMode;
  
  bool normalize;
  KernelType spectrumKernel;

  char * weights_file;
  char * init_vit_file;
  
  int optimisation_method;
  
  int init_weights;
  
};

/**
 * Show the program's usage.
 */

void showUsage ()
{
  cerr << "Usage : moleculesWeightedRegression trainset_file testset_file [options]" << endl;
  cerr << "options:" << endl;
	
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;
  cerr << "\t5 -- Complete Gaussian Kernel" << endl;

  cerr << "-l lambda : Set the parameter of the KRR" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel " << endl;
  cerr << "-m mu : Set the paramater of the W regularisation term" << endl;
  cerr << "-o file : Export des poids" << endl;
  cerr << "-t int : Optimisation Method " << endl;
  cerr << "-v int : Initialisation des poids " << endl;
  cerr << "\t 1 : Acyclic VIT par Forward Selection" << endl;
  cerr << "\t -1 : Fichier de VIT" << endl;
  cerr << "-i file : Fichier Initialisation des VIT" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-N : Normalize Treelet Kernel" << endl;
}

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
  options->lambda = 1.0;
  options->mu = 1.0;
  options->trainset_file = argv[1];
  options->testset_file = argv[2];
  options->quietMode = false;
  options->normalize = false;
  options->spectrumKernel = InnerProductKernelType;
  options->weights_file = (char*)"";
  options->init_vit_file = (char*)"";
  options->optimisation_method = 0;
  options->init_weights = 0;
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
      else if (strncmp(argv[i], "-m", 2) == 0)
	{
	  options->mu = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-N", 2) == 0)
	{
	  options->normalize= true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->spectrumKernel = (KernelType)(atoi(argv[i+1]));
	  i+=2;
	}
      else if (strncmp(argv[i], "-o", 2) == 0)
	{
	  options->weights_file = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->optimisation_method = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-v", 2) == 0)
	{
	  options->init_weights = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-i", 2) == 0)
	{
	  options->init_vit_file = argv[i+1];
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
	
  TreeletKernel* kgraph = NULL;
  MoleculeGraph::initTable();
  
  char * trainset_file_copy = new char[strlen(options.trainset_file)];
  strcpy(trainset_file_copy,options.trainset_file);
   
  MoleculesDataset trainset (options.trainset_file);
  MoleculesDataset testset (options.testset_file);
  MoleculesDataset dataset (trainset_file_copy);
    
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
    case CompleteGaussianKernelType:
      cout << "Complete Gaussian Kernel";
      break;
    default:
      //devrait pas arriver
      break;
    }
  cout << endl;
  
  trainset.computeLabeledSpectrums();
  testset.computeLabeledSpectrums();

  if(options.normalize)
    {
      trainset.normalizeLabeledSpectrums();
      testset.normalizeLabeledSpectrums();
    }
  
  kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
  if(!options.quietMode){
    cout << "Dataset : " << options.trainset_file << endl;
    cout << "Number of molecules : " << trainset.size() << endl;
    cout << endl;
  }
  
  // We compute the gram matrix of the training set
  // trainset.computeGramMatrix(kgraph, true);
  // testset.computeGramMatrix(kgraph, true);

  // if(!options.quietMode)
  //   cout << "Gram Matrix Positive definite ? ! " << dataset.isGramMatrixPD() << endl;;
  
  // cimg_library::CImg<double> gram = trainset.getGramMatrix ();
  
  // cout << "Affichage des similaires :"<< endl;
  // for(int i=0;i<dataset.size();i++)
  //   for(int j=0;j<dataset.size();j++)
  //     if (gram(i,j) == 1 && (i!=j))
  // 	cout << basename(dataset.getMoleculeFilename(i)) << "(" << i << ")" << " = " 
  // 	     << basename(dataset.getMoleculeFilename(j)) << "(" << j << ")" << endl;
  

  
  WeightedRegression optim(&trainset);
  optim._epsilon = optim._epsilon_alpha = 0.001;
  optim._epsilon_w = 0.00001;

  optim._sigma = options.sigma;
  optim._lambda = options.lambda;
  optim._mu = options.mu;

  optim._normalize = false;
  CImg<double> w0;
  switch(options.init_weights)
    {
    case 1:
      w0 = optim.vit_to_weights(getAcyclicVIT_2());
      break;
    case -1:
      w0 = optim.vit_to_weights(getAcyclicFromFile(options.init_vit_file));
      break;
    default:
      //pas de w0
      ;
      
    }
  
  if(options.optimisation_method == 0)
    optim.computeMinimisation(kgraph, w0);
  else if (options.optimisation_method == 1)
    optim.computeMinimisationV2(kgraph, w0);
  else
    {
      cerr << "Unknown Optimisation Method"  << endl;
      showUsage();
      exit(1);
    }
  
  //Prediction
  CImg<double> weights = optim.getOptimalWeights();
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  //Conversion Cimg<->Map
  map<string, double,bool (*)(string, string)> weights_map[SIZE_SPECTRUM];
  for(int n=0;n<SIZE_SPECTRUM;n++)
    {
      weights_map[n] = map<string, double,bool (*)(string, string)> (string_utils::keyStringComp);
      //vit_list[0] = vector<string> ();
    }
  
  for(unsigned int k=0;k<weights.size();k++)
    {
      weights_map[optim._k_to_treelet[k].second].
      	insert(pair<string,double>(optim._k_to_treelet[k].first,weights[k]));
      
      vit_list[optim._k_to_treelet[k].second].push_back(optim._k_to_treelet[k].first);
    }

  kgraph->selectTreelets(vit_list);
  kgraph->weightTreelets(weights_map);
  trainset.computeGramMatrix(kgraph, false);
  testset.computeGramMatrix(kgraph, false);

  // CImg<double> alpha = optim.getOptimalAlpha();
  for(int i=0;i<weights.height();i ++)
    cout << weights(0,i) << endl;
  // int N = testset.size();
  /**Test**/
  dataset.computeLabeledSpectrums();
  TreeletKernel * kgraph2 = new TreeletKernel(options.spectrumKernel,options.sigma);
  kgraph2->selectTreelets(vit_list);
  kgraph2->weightTreelets(vit_to_weight(vit_list));
  //kgraph2->weightTreelets(weights_map);
  dataset.computeGramMatrix(kgraph2, false);
  KRegression* kr = new KernelRidge(kgraph2, &dataset, options.lambda);   
  /**/ 
  double erreur_without_weigths = 0.0;
  double erreur_with_weigths = 0.0;
  for (unsigned int i=0; i<testset.size(); ++i)
    {
      /*Version sans poids*/
      Collection* col = testset.getCollection(i);
      double y_estim = (*kr) (col);
      cout << "estim " << i << ": " << y_estim << ", " << testset.getParameter(i) << endl;
      erreur_without_weigths += pow(y_estim - testset.getParameter(i),2);
      
      
      /*Version avec poids*/
      CImg<double> K (1,trainset.size());
      CImg<double> y (trainset.size(),1);
      
      for(unsigned int n=0;n<trainset.size();n++)
      	{
	      double k = (*kgraph)(col,trainset.getCollection(n));
	      double a = (*kgraph) (col, col);
	      double b = (*kgraph) (trainset.getCollection(n),trainset.getCollection(n));
	      k /= sqrt(a*b);	
	      
	      K(0,n) = k;
	      y(n,0) = trainset.getParameter(n);
      	}
      CImg<double> tmp = trainset.getGramMatrix(true) + 
      	CImg<double>(trainset.size(),trainset.size()).identity_matrix()*options.lambda;
      CImg<double> tmp2 = tmp.invert();
      CImg<double> res = y*tmp2*K;
      double y_estim2 = res(0,0);
      erreur_with_weigths += pow(y_estim2 - testset.getParameter(i),2);      
      cout << "estim KRR " << i << ": " << y_estim2 << ", " << testset.getParameter(i) << endl;
    }

  cout << "EQM sans pondération : " << erreur_without_weigths / testset.size() << endl;
  cout << "EQM avec pondération : " << erreur_with_weigths / testset.size() << endl;
  
  if(strlen(options.weights_file) != 0)
    {
      ofstream f (options.weights_file);
      for(unsigned int k=0;k<weights.size();k++)
	{
	  double weight = weights[k];
	  int treelet_type = optim._k_to_treelet[k].second;
	  string code_treelet = optim._k_to_treelet[k].first;
	  f << treelet_type << "-" << code_treelet << " " << weight << endl;
	  }
      f.close();
    }

  trainset.destroy();
  testset.destroy();	
  delete kgraph;
  return 0;
}

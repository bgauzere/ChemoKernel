/*
 * @file moleculesRegression10percent.cpp
 *
 * @description This program performs a Kernel Ridge regression on molecules
 *              datasets using 90% of the dataset as traing set end the remaining 10% as the test set.
 * 
 * @author Benoit Gauzere <benoit.gauzere@ensicaen.fr>
 *
 * @version 1.0.0 (2010-10-26)
 */

#include <iostream>
#include <pandore.h>
#include <cfloat>
#include "MoleculeGraph.h"
#include "ShapeBagTrails.h" 
#include "MoleculesDataset.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"

#include "GraphKernel.h"
#include "TrailKernel.h"
#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "KRegression.h"
#include "GraphletCountKernel.h"
#include "KernelRidge.h"
#include "LaplacianKernelOriginal.h"
#include "KEditDistance.h"
#include "Kashima.h"
#include "TreeletKernel.h"
#include "utils.h"
#include "string_utils.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"


using namespace std;
using namespace pandore;
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

bool isVIT(string treelet, vector<string> vit_list){
  for(unsigned int i = 0;i< vit_list.size();i++)
    if(! vit_list[i].compare(treelet))
      return true;
  return false;
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
vector<string> * getStepWiseVIT (MoleculesDataset & dataset,  TreeletKernel* kgraph, double lambda)
{	
  dataset.computeLabeledSpectrums();
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  int N = dataset.size();
  double std_dev = DBL_MAX;
  double std_dev_current = 80;
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  // for(int i=0;i<SIZE_SPECTRUM;i++)
  //   vit_list[i] = vector<string>;
  vector<string> * candidate_vit_list = new vector<string>[SIZE_SPECTRUM];
  //XXX: Vieux Hack de sioux pour eviter le produit scalaire de deux trucs vides 
  
  cout << (int)distribution[0].begin()->first[0]<< endl;
  candidate_vit_list[0].push_back(distribution[0].begin()->first);  
  
  int nb_regressions = 0;
  int nb_patterns = 0;
  int nb_vit = 0;
  //XXX: Vérif de la complexité
  //while(std_dev_current < std_dev)
  for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      treelet_spectrum::iterator it = distribution[treelet_type].begin();
      for(;it != distribution[treelet_type].end();it ++)
	{
	  nb_patterns ++;
	}
    }
  while(std_dev_current < std_dev)
    //while(nb_vit < nb_patterns)
    { //On s'arrete quand on progresse plus, (évolution a plotter ds un 2nd temps)
      nb_vit ++;
      std_dev = std_dev_current;
      std_dev_current = DBL_MAX;
      for(int i=0;i<SIZE_SPECTRUM;i++)
	vit_list[i] = vector<string>(candidate_vit_list[i]);
      //vit_list = candidate_vit_list;
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	{
	  treelet_spectrum::iterator it = distribution[j].begin();
	  for(;it != distribution[j].end();it ++)
	    //Suppression des peu fréquents
	    //if(it->second > LIMIT_LOW_FREQUENCY)
	    if(! isVIT(it->first,vit_list[j]))
	      {//Regression avec chaque treelet != vit_list
		//Creation de la vit_list candidate
		vector<string> * tmp_vit_list = new vector<string>[SIZE_SPECTRUM];
		for(int i=0;i < SIZE_SPECTRUM;i++)
		  tmp_vit_list[i] = vector<string>(vit_list[i]);
		tmp_vit_list[j].push_back(it->first);
		kgraph->selectTreelets(tmp_vit_list);
		kgraph->weightTreelets(vit_to_weight(tmp_vit_list));
		dataset.computeGramMatrix(kgraph, false);
		//Regression
		KRegression* kr = new KernelRidge (kgraph, &dataset, lambda);
		double std_dev_tmp;
		double err_quad = 0.0;
		for (int i=0; i<N; ++i)
		  {
		    double bp_exp = dataset.getParameter(0);
		    Collection* col = dataset.getCollection(0);
		    dataset.delete_first();
		    double bp = (*kr) (col);
		    double err = bp_exp - bp;
		    err_quad += err * err;
		    dataset.add (col, bp_exp, kgraph);
		  }
		std_dev_tmp = sqrt(err_quad / N);
		nb_regressions ++;
		//cout << it->first << " : " << std_dev_tmp << endl;
		if(std_dev_tmp < std_dev_current)
		  {
		    std_dev_current = std_dev_tmp;
		    for(int i=0;i<SIZE_SPECTRUM;i++)
		      candidate_vit_list[i] = vector<string>(tmp_vit_list[i]);
		  }
		else
		  {
		    
		  }
		delete [] tmp_vit_list;
	      }
	}
      cout << nb_vit << " " << std_dev_current << endl;
      //Ici, on a un nouveau vit + un nouveau ecart type 
    }
  return  vit_list;
}

vector<string> * getBackwardStepWiseVIT (MoleculesDataset & dataset,  TreeletKernel* kgraph, double lambda)
{	
  cout << "Initialisation"  << endl;
  //  dataset.computeLabeledSpectrums();
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  int N = dataset.size();
  double std_dev = DBL_MAX;
  double std_dev_current = 80;
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  vector<string> * candidate_vit_list = new vector<string>[SIZE_SPECTRUM];
  int nb_patterns = 0;
  //Initialisation de la vit list avec tous les treelets
  for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      treelet_spectrum::iterator it = distribution[treelet_type].begin();
      for(;it != distribution[treelet_type].end();it ++)
	{
	  candidate_vit_list[treelet_type].push_back(it->first);
	  nb_patterns ++;
	}
    }
  
  
  int nb_vit = nb_patterns;
  while(std_dev_current < std_dev) /*On s'arrete quand on progresse plus, (évolution a plotter ds un 2nd temps)*/
    //(nb_vit != 0)
    //while(nb_vit != 56) // >0
    { 
      //cout << "step " << nb_patterns - nb_vit << endl;
      int nb_regressions = 0;
      std_dev = std_dev_current;
      std_dev_current = DBL_MAX;
      for(int i=0;i<SIZE_SPECTRUM;i++)
	vit_list[i] = vector<string>(candidate_vit_list[i]);
      
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	{
	  treelet_spectrum::iterator it = distribution[j].begin();
	  for(unsigned int i=0;i<candidate_vit_list[j].size();i++) 
	    {
	      //On supprime candidate_vit_list[j][i]
	      //Creation de la vit_list candidate
	      vector<string> * tmp_vit_list = new vector<string>[SIZE_SPECTRUM];
	      for(int treelet_type=0;treelet_type < SIZE_SPECTRUM;treelet_type++)
		tmp_vit_list[treelet_type] = vector<string>(vit_list[treelet_type]);
	      
	      tmp_vit_list[j].erase(tmp_vit_list[j].begin()+i);
	      kgraph->selectTreelets(tmp_vit_list);
	      kgraph->weightTreelets(vit_to_weight(tmp_vit_list));
	      
	      dataset.computeGramMatrix(kgraph, false);
	      //Regression
	      KRegression* kr = new KernelRidge (kgraph, &dataset, lambda);
	      double std_dev_tmp;
	      double err_quad = 0.0;
	      for (int i=0; i<N; ++i)
		{
		  double bp_exp = dataset.getParameter(0);
		  Collection* col = dataset.getCollection(0);
		  dataset.delete_first();
		  double bp = (*kr) (col);
		  double err = bp_exp - bp;
		  err_quad += err * err;
		  dataset.add (col, bp_exp, kgraph);
		}
	      std_dev_tmp = sqrt(err_quad / N);
	      nb_regressions ++;
	      if(std_dev_tmp < std_dev_current)
		{
		  std_dev_current = std_dev_tmp;
		  for(int treelet_type=0;treelet_type<SIZE_SPECTRUM;treelet_type++)
		    candidate_vit_list[treelet_type] = vector<string>(tmp_vit_list[treelet_type]);
		}
	      else
		{
		  
		}
	      delete [] tmp_vit_list;
	    }
	  
	}
      // for(int treelet_type=0;treelet_type<SIZE_SPECTRUM;treelet_type++)
      // 	for(int i=0;i<candidate_vit_list[treelet_type].size();i++)
      // 	  cout << treelet_type << candidate_vit_list[treelet_type][i] << endl;
      cout << endl << nb_vit << " " << std_dev_current << endl;
      //Ici, on a un nouvel ensemble de vit + un nouveau ecart type 
      nb_vit --;
    }
  return  vit_list;
}


double predict (MoleculesDataset & trainset, GraphKernel* kgraph, double lambda, Collection* col)
{	
  KRegression* kr = new KernelRidge (kgraph, &trainset, lambda);
  return (*kr) (col);
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
  double mu;
  TreeletKernel::KernelType spectrumKernel;
  bool normalize;
  int top;
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
  cerr << "\t4 -- Graphlet Count Kernel" << endl;
  cerr << "\t5 -- Random Walks Kernel" << endl;
  cerr << "\t6 -- Treelet Kernel" << endl;

  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;
  cerr << "\t5 -- Complete Gaussian Kernel" << endl;
	
  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-N : Normalize Treelet Kernel" << endl;
  cerr << "-t t: Treelet selection" << endl;
  cerr << "\t t = 1 : Forward Selection" << endl;
  cerr << "\t t = 2 : Backward Selection" << endl;
}

void readOptions (int argc, char** argv, Options* options)
{
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  if (argc < 3)
    {
      usage(argv[0]);
      exit(1);
    }
	
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
  options->mu = 1.0;
  options->spectrumKernel = TreeletKernel::InnerProductKernelType;
  options->normalize = false;
  options->top = 0;
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
      else if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->spectrumKernel = (TreeletKernel::KernelType)(atoi(argv[i+1]));
	  i+=2;
	}
      else if (strncmp(argv[i], "-N", 2) == 0)
	{
	  options->normalize= true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->top = atoi(argv[i+1]);
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

class TestSetSplit {

public:
  
  Options options;
  double * bp_pred;
  int step ; // 1/10  molecule
  const MoleculesDataset d;
  
  TestSetSplit(Options options,  double * bp_pred, const MoleculesDataset& d) : 
    options(options), bp_pred(bp_pred), d(d){
    step = 10;
        
  }

  void operator() ( const tbb::blocked_range<int>& r ) const {
    
    for (int k = r.begin(); k != r.end(); ++k) {
      cout << k << endl;
      ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
      GraphEditDistance* edit = new MoleculeGraphEditDistance;
      GraphKernel* kgraph = NULL;
      MoleculesDataset * dataset = new MoleculesDataset (d);
      int N = dataset->size();
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
      switch(options.kernel)
	{
	case 0:
	  kgraph = new KMean(sbt, NULL);
	  break;
	case 1:
	  kgraph = new KWMean (sbt, NULL, options.sigma);
	  break;
	case 2:
	  kgraph = new LaplacianKernelOriginal (edit, *dataset, options.sigma, options.regularization, options.lambda);
	  break;
	case 3:
	  kgraph = new KEditDistance (edit, options.sigma);
	  break;	
	case 4:  
	  kgraph = new GraphletCountKernel(options.sigma,options.mu);
	  coeffs = new double[SIZE_SPECTRUM];
	  for(int i = 0;i < SIZE_SPECTRUM;i++)
	    coeffs[i] = dataset->getCorrelationCoeff(i);
	  ((GraphletCountKernel*)kgraph)->setCoeffs(coeffs,SIZE_SPECTRUM);
	  break;
	case 5:  
	  kgraph = new RandomWalkKernel(new Kashima(),options.lambda,1.0);
	  break;			
	case 6:
	  kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
	  if(options.top == 0){
	    vector<string> * vit_list = getAllAsVIT(*dataset);
	    ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	    ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	  }else if (options.top == -1){//Forward Selection
	    MoleculeGraph util;
	    vector<string> * vit_list = getStepWiseVIT(*dataset,(TreeletKernel*)kgraph,options.alpha);
	    cout << "VIT List : " << endl;
	    for(int i=0;i<SIZE_SPECTRUM;i++)
	      for(int j=0;j<vit_list[i].size();j++)
		cout << "G" << i << "-" << vit_list[i][j] << endl;
	    ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	    ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	  }else if (options.top == -2){//Backward Elimination
	    //MoleculeGraph util;
	    vector<string> * vit_list = getBackwardStepWiseVIT(*dataset,(TreeletKernel*)kgraph,options.alpha);
	    cout << "VIT List : " << endl;
	    for(int i=0;i<SIZE_SPECTRUM;i++)
	      for(unsigned int j=0;j<vit_list[i].size();j++)
		cout << "G" << i << "-" << vit_list[i][j] << endl;
	    ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	    ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	  }
	  break;
	default:
	  cerr << "Error : Invalid kernel" << endl;
	  exit(1);
	}
      dataset->computeGramMatrix(kgraph, false);
    
      //Test molecules prediction
      list<int>::iterator it = to_test.begin();
      list<Collection*>::iterator it_col = to_test_Collections.begin();    
      list<double>::iterator it_param = to_test_Parameters.begin();
      
      for(;it != to_test.end();++it, ++it_col, ++it_param)
	bp_pred[*it] = predict(*dataset, kgraph, options.alpha,*it_col);
      
      //Dataset cleaning
      // dataset->destroy();
      // delete kgraph;
    }
  }
};


int main (int argc, char** argv)
{	
  //Display command for archiving results
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  
  //Initializations
  Options options;
  readOptions(argc, argv, &options);
  
  MoleculeGraph::initTable();
  MoleculesDataset dataset(options.trainset_path, options.dataset_file);
  int N = dataset.size();
  int step = 10; // 1/10  molecule
  //Stats data
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;
  double * bp_exps = new double[N];
  double * bp_pred = new double[N];
  //Pre compute Graphlet Count Kernel
  if(options.kernel == 4){
    dataset.computeSpectrums();
    dataset.computeGraphletCorrelation();
  }  
  //Pre compute Labeled Treelet distribution
  if(options.kernel == 6){
    dataset.computeLabeledSpectrums();
    if(options.normalize)
      dataset.normalizeLabeledSpectrums();
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, step),TestSetSplit(options,bp_pred,dataset));

  for(int i = 0; i < N; i ++){
    bp_exps[i] = dataset.getParameter(i);
  }
  
  err_quad =0;
  err_moy =0;
  bp_mean =0;

  for(int i=0;i<N;i++){
    cout << "Molecule " << i << " : " << bp_pred[i] << " (" << bp_exps[i] << ")" << endl;
    double err =  bp_exps[i] - bp_pred[i];
    err_quad +=  err * err;
    err_moy += abs(err);
    bp_mean += bp_exps[i]; 
  }
  
  bp_mean /= N;
  std_dev = sqrt(err_quad / N); 
  err_moy /= N;

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

  return 0;
}


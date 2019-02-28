/*
 * @file moleculesRegressionTrainAndTestSet.cpp
 *
 * @description This program performs a Kernel Ridge regression or a SVR on two molecules
 *              datasets using the first dataset as traing set end the other as the test set.
 * 
 * @author Benoit Gauzere <benoit.gauzere@ensicaen.fr>
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *
 * @version 1.0.0 (2012-11-08)
 */




#include <iostream>
#include <pandore.h>
#include <cfloat>
#include <sys/times.h>
#include <sys/time.h>
#include <fstream>
#include <cassert>
#include <sstream>

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
#include "SVM.h"
#include "LaplacianKernelOriginal.h"
#include "KEditDistance.h"
#include "Kashima.h"
#include "TreeletKernel.h"
#include "KPrototypeEditDistance.h"
#include "utils.h"
#include "string_utils.h"
#include "TreeletCycleKernel.h"
#include "ContractedCycleKernel.h"
#include "InterTreeletKernel.h"
#include "GaussianKernel.h"
#include "InnerProductKernel.h"
#include "IntersectionKernel.h"
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
bool always_true(double i)
{
  return true;
}




vector<string> * getSizeTreeletsasVIT (MoleculesDataset & dataset, int taille){
  // initializer_list<int> mylist;
  // int nb_treelets = 0;
  // switch (taille){
  // case 1:
  //   mylist = {0};
  //   nb_treelets = 1;    
  //   break;
  // case 2:
  //   mylist = {0,1};
  //   nb_treelets = 2;    
  //   break;
  // case 3:
  //   mylist = {0,1,2};
  //   nb_treelets = 3;    
  //   break;
  // case 4:
  //   mylist = {0,1,2,3,6};
  //   nb_treelets = 5;    
  //   break;
  // case 5:
  //   mylist = {0,1,2,3,4,6,7,8};
  //   nb_treelets = 8;    
  //   break;
  // case 7:
  //   mylist = {0,1,2,3,4,5};
  //   nb_treelets = 6;    
  //   break;
  // default:
  //   mylist = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
  //   nb_treelets = 14;    
  //   ;
  // }
  // vector<int> list_struct(mylist.begin(),mylist.end());
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


vector<string> * getSizeContractedCycleTreeletsasVIT(MoleculesDataset & dataset, int taille){
  // initializer_list<int> mylist;
  // int nb_treelets = 0;
  // switch (taille){
  // case 1:
  //   mylist = {0};
  //   nb_treelets = 1;    
  //   break;
  // case 2:
  //   mylist = {0,1};
  //   nb_treelets = 2;    
  //   break;
  // case 3:
  //   mylist = {0,1,2};
  //   nb_treelets = 3;    
  //   break;
  // case 4:
  //   mylist = {0,1,2,3,6};
  //   nb_treelets = 5;    
  //   break;
  // case 5:
  //   mylist = {0,1,2,3,4,6,7,8};
  //   nb_treelets = 8;    
  //   break;
  // case 7:
  //   mylist = {0,1,2,3,4,5};
  //   nb_treelets = 6;    
  //   break;
  // default:
  //   mylist = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
  //   nb_treelets = 14;    
  //   ;
  // }
  // vector<int> list_struct(mylist.begin(),mylist.end());
  // treelet_spectrum * distribution = dataset.getContractedCycleTreeletDistribution(&always_true);
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


vector<string> * getAzoteAsVIT (MoleculesDataset & dataset){
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
	  string label = it->first;
	  //On ne garde que les azotes
	  if(label.find(char(7)) != string::npos){
	    vit_list[treelet_type].push_back(it->first);
	    nb_patterns ++;
	  }
	}
    }
  return vit_list;
}

vector<string> * getCycleAzoteAsVIT (MoleculesDataset & dataset){
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
	  string label = it->first;
	  //On ne garde que les azotes
	  if(label.find(char(7)) != string::npos){
	    vit_list[treelet_type].push_back(it->first);
	    nb_patterns ++;
	  }
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
  dataset.computeLabeledSpectrums();
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


double predict (MoleculesDataset & trainset, GraphKernel* kgraph, double lambda, Collection* col,int methode,double c, double epsilon)
{ 
  double res;
  if(methode==2)
    {
      SVM svm (trainset, kgraph, SVM::EPSILON_SVR,c,0.5,epsilon);
      res=svm.predict(col);
    }
  else
    {
      KRegression* kr = new KernelRidge (kgraph, &trainset, lambda);
      res=(*kr) (col);
    }
 
  return res;
}

vector<string> * getContractedCycleAzoteAsVIT (MoleculesDataset & dataset){
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
	  string label = it->first;
	  //On ne garde que les azotes
	  if(label.find(char(7)) != string::npos){
	    vit_list[treelet_type].push_back(it->first);
	    nb_patterns ++;
	  }
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


struct Options
{
  char* dataset_path;
  char* testset_file;
  char* trainset_file;
  int kernel;
  int regMethode;
  double sigma;
  int regularization;
  double c;
  double lambda;
  double alpha;
  double keptTrails;
  double trailsLength;
  int editdistance;
  bool quietMode;
  double mu;
  KernelType spectrumKernel;
  KernelType spectrumKernelCycle;
  KernelType spectrumKernelContractedCycle;
  int nb_neighbours;
  int step;
  bool normalize;
  int top;
  char * weights_file;
  KPrototypeEditDistance::PrototypeSelectionType prototypeSelection;
  double cs,ci; //costs for substitution and insertion/deletion for edit distance
  double sigma_inter;
  int foldCV;
  int labelCycles;
  double svm_epsilon;
  int size_max;
};

void usage (char * s)
{
  cerr << "Usage : " << s << " dataset_path trainset_file testset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-k kernel_type : Set the type of kernel function (default 0)" << endl;
  cerr << "\t0 -- Kernel Mean" << endl;
  cerr << "\t1 -- Kernel Weighted Mean" << endl;
  cerr << "\t2 -- Graph Laplacian Kernel" << endl;
  cerr << "\t3 -- Gaussian Kernel" << endl;
  cerr << "\t4 -- Graphlet Count Kernel" << endl;
  cerr << "\t5 -- Random Walks Kernel" << endl;
  cerr << "\t6 -- Treelet Kernel" << endl;
  cerr << "\t10 -- Prototype Edit Distance Kernel" << endl;
  cerr << "\t11 -- Cycle Kernel" << endl;
  cerr << "\t13 -- Inter Treelet Kernel - Approximate Edit Distance" << endl;
  cerr << "\t14 -- Inter Treelet Kernel - Laplacian Kernel" << endl;
  cerr << "\t15 -- Inter Treelet Kernel - Exact Connexe Edit Distance (Regul)" << endl;
  cerr << "\t19 -- Treelet Kernel on Contracted Cycle Hypergraph" << endl;
  
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;
  cerr << "\t5 -- Complete Gaussian Kernel" << endl;

  cerr << "-P method : Set the method to select the prototypes" << endl;
  cerr << "-r nb proto : using -k 10, set the number of prototypes" << endl;

  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-eps epsilon : Set the parameter for the size of epsilon tube of SVM Regression (default 0.1)" << endl;
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
  cerr << "\t t = 7 : Selection des Azotes" << endl;
  cerr << "-ci : Node/Edge Insertion/Deletion cost" << endl;
  cerr << "-cs : Node/Edge substitution cost" << endl;
  cerr << "-sigma_inter : Sigma pour la similarité inter treelets" << endl;

  cerr << "-reg regression : Choisis différentes méthodes de régression" << endl;
  cerr << "\t1 -- Kernel Ridge" << endl;
  cerr << "\t2 -- SVM" << endl;

  cerr << "-c C SVM : Paramètre C du SVM" << endl; 
  cerr << "\t1 -- Parcours par couche" << endl;
  cerr << "\t2 -- Parcours par branche" << endl;

  cerr << "-fo Fold Value use for the CV (the CV is done to calcul q) if set to 0 CV is done by LOO" << endl;

  cerr << "-lC Modifications of label"<< endl;
  cerr << "\t0 -- No modifications" << endl;
  cerr << "\t1 -- Number of relevant cycles add to the label of an atom" << endl;
  cerr << "\t2 -- Size of each relevant cycles add to the label of an atom" << endl;
  cerr << "\t-size -- Maximum size of considered treelets" << endl;
}

void readOptions (int argc, char** argv, Options* options)
{
  if (argc < 4)
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
  options->dataset_path = argv[1];
  options->trainset_file = argv[2];
  options->testset_file = argv[3];
  options->mu = 1.0;
  options->quietMode = false;
  options->spectrumKernelCycle = InnerProductKernelType;
  options->spectrumKernel = InnerProductKernelType;
  options->spectrumKernelContractedCycle = InnerProductKernelType;  
  options->normalize = false;
  options->top = 0;
  options->nb_neighbours = 10;
  options->step = 3;
  options->weights_file = NULL;
  options->ci = 7;
  options->cs = 1;
  options->sigma_inter = 1;
  options->regMethode=1;
  options->c=1.0;
  options->foldCV=0;
  options->labelCycles=0;
  options->svm_epsilon=0.1;
  options->size_max=6;

  int i=4;
  while (i<argc)
    {
      if (strncmp(argv[i], "-lC", 3) == 0)
	{
	  options->labelCycles = atoi(argv[i+1]);
	  i+=2;
	} 
      else if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->kernel = atoi(argv[i+1]);
	  i+=2;
	} 
      else if (strncmp(argv[i], "-reg", 4) == 0)
	{
	  options->regMethode = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-fo", 3) == 0)
	{
	  options->foldCV = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strcmp(argv[i], "-size") == 0)
	{
	  options->size_max = atoi(argv[i+1]);
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
      else if (strcmp(argv[i], "-eps") == 0)
	{
	  options->svm_epsilon = atof(argv[i+1]);
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
      else if (strncmp(argv[i], "-b", 2) == 0)
	{
	  options->nb_neighbours = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-x", 2) == 0)
	{
	  options->step = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-w", 2) == 0)
	{
	  options->weights_file = argv[i+1];
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
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atoi(argv[i+1]);
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
  //Display command for archiving results
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  
  //Initializations
  Options options;
  readOptions(argc, argv, &options);
  ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
  GraphKernel* kgraph = NULL;
  GraphEditDistance* edit = new MoleculeGraphEditDistance;
  
  MoleculeGraph::initTable();
  MoleculesDataset * dataset = new MoleculesDataset (options.dataset_path, options.trainset_file);
  //Pre compute Graphlet Count Kernel
  if(options.kernel == 4){
    dataset->computeSpectrums();
    dataset->computeGraphletCorrelation();
  }  
  //Pre compute Labeled Treelet distribution
  if(options.kernel >= 6 && options.kernel != 16){
    dataset->computeLabeledSpectrums();
    if(options.normalize)
      dataset->normalizeLabeledSpectrums();
  }

  if(options.kernel == 19){
    dataset->computeLabeledSpectrums();
    dataset->computeCCHypergraph();
    dataset->computeContractedCycleSpectrums();
  }

  int N = dataset->size();
  
  long mean_prediction_time = 0;
  
  MoleculesDataset * testset = new MoleculesDataset (options.dataset_path, options.testset_file);


  if(options.kernel == 4){
    testset->computeSpectrums();
    testset->computeGraphletCorrelation();
  }  

  if(options.kernel >= 6 && options.kernel != 16){
    testset->computeLabeledSpectrums();
    if(options.normalize)
      testset->normalizeLabeledSpectrums();
  }
  
  if(options.kernel == 19){
    testset->computeLabeledSpectrums();
    testset->computeCCHypergraph();
    testset->computeContractedCycleSpectrums();
  }

  int N2 = testset->size();

  //Stats data
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;
  double bp_exp_mean = 0.0;
  double * bp_exps = new double[N2];
  double * bps = new double[N2];
  int nb_molecules = 0;
  for(int i = 0; i < N2; i ++){
    bp_exps[i] = testset->getParameter(i);
    bp_exp_mean += bp_exps[i];
  }
  bp_exp_mean /= N2;

  list<int> to_test;
  vector<Collection *> to_test_Collections;
  list<double> to_test_Parameters;
  vector<MoleculeGraph*> to_test_Molecules;
  for(int i = 0; i < N2; i ++){
    to_test.push_back(i);
    to_test_Collections.push_back(testset->getCollection(i));
    to_test_Parameters.push_back(testset->getParameter(i));
    to_test_Molecules.push_back(testset->getMoleculeGraph(i));
  }
    
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

      if(options.top == 0){
	vector<string> * vit_list = getSizeTreeletsasVIT(*dataset, options.size_max);
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
      }else if (options.top == -1){//Forward Selection
	MoleculeGraph util;
	vector<string> * vit_list = getStepWiseVIT(*dataset,(TreeletKernel*)kgraph,options.alpha);
	cout << "VIT List : " << endl;
	for(int i=0;i<SIZE_SPECTRUM;i++)
	  for(unsigned int j=0;j<vit_list[i].size();j++)
	    cout << "G" << i << "-" << vit_list[i][j] << endl;
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
      }else if (options.top == -2){//Backward Elimination
	//MoleculeGraph util;
	vector<string> * vit_list = getBackwardStepWiseVIT(*dataset,(TreeletKernel*)kgraph,options.alpha);
	cout << "VIT List : " << endl;
	for(int i=0;i<SIZE_SPECTRUM;i++)
	  for(int j=0;j<vit_list[i].size();j++)
	    cout << "G" << i << "-" << vit_list[i][j] << endl;
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
      }
      else if (options.top == 7){//Selection des azotes
	//MoleculeGraph util;
	cout << "Selection des treelets d'azote" << endl;
	vector<string> * vit_list = getAzoteAsVIT(*dataset);
	cout << "VIT List : " << endl;
	for(int i=0;i<SIZE_SPECTRUM;i++)
	  for(int j=0;j<vit_list[i].size();j++)
	    cout << "G" << i << "-" << MoleculeGraph::translateTreeletCode(vit_list[i][j]) << endl;
	((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
      }
      

      break;
    case 10:
      cout << "Kernel : Prototype Edit Distance" << endl;
      kgraph = new KPrototypeEditDistance(*dataset, options.regularization,options.prototypeSelection);
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
      dataset->computeRelevantCycles();
      dataset->computeCyclesSpectrums();
      
      kgraph = new CycleKernel(options.spectrumKernelCycle,options.sigma);
      
      if(options.weights_file == NULL){
	((CycleKernel*)kgraph)->selectTreelets(getCyclesAllAsVIT(*dataset));
	((CycleKernel*)kgraph)->weightTreelets(vit_to_weight(getCyclesAllAsVIT(*dataset)));
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
      dataset->getCycleTreeletDistribution(&always_true);
      cout << "Nb Cycle Treelets : " << dataset->getNbCyclesTreelets() << endl;
      break;
    case 13:      
      {
	dataset->computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset->getTreeletDistribution(&always_true);
	int nb_treelets =  dataset->getNbTreelets();
	dataset->computeSpecialVectors(&distrib, nb_treelets);
	for(int i = 0;i<to_test_Molecules.size();i++)
	  to_test_Molecules[i]->computeSpecialLabeledTreeletVector(&distrib,nb_treelets);
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
	GraphEditDistance * ed_treelets = new MoleculeGraphEditDistanceMCS(options.cs, options.ci);
	GraphKernel * kinter_treelets = new KEditDistance (ed_treelets, options.sigma_inter);

	kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
					&distrib,
					nb_treelets);
      }
      break;
    case 14:      
      {
	dataset->computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset->getTreeletDistribution(&always_true);
	int nb_treelets =  dataset->getNbTreelets();
	dataset->computeSpecialVectors(&distrib, nb_treelets);
	for(int i = 0;i<to_test.size();i++)
	  to_test_Molecules[i]->computeSpecialLabeledTreeletVector(&distrib,nb_treelets);
	cout << "Special Vectors computed " << endl;
	
	Dataset dataset_treelets;
	Collection ** col_treelets = new Collection* [nb_treelets];
	int cur_treelet = 0;
	for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
	  {
	    treelet_spectrum::iterator it = distrib[i].begin();
	    for(;it != distrib[i].end();it++,cur_treelet++){
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
	cout << "InterTreeletKernel Initialized" << endl;
      }
      break;
    case 15:      
      {
	dataset->computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset->getTreeletDistribution(&always_true);
	int nb_treelets =  dataset->getNbTreelets();
	cout << nb_treelets<< endl;
	dataset->computeSpecialVectors(&distrib, nb_treelets);
	cout << "Special Vectors computed " << endl;
	for(int i = 0;i<to_test.size();i++)
	  to_test_Molecules[i]->computeSpecialLabeledTreeletVector(&distrib,nb_treelets);

	
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
	  

	kgraph = new InterTreeletKernel(NULL,kinter_graphs,
					&distrib,
					nb_treelets);
	
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
	    }
	  }
        CImg<double> sim_treelets(nb_treelets,nb_treelets,1,1,0.0);
	// for (unsigned int i=0; i<nb_treelets; ++i)
	//   for(int j=0;j<nb_treelets;j++)
	//     sim_treelets(i,j) = 1;
	
	TreeletEditDistance ed(options.cs,options.ci);	
	TreeletEditDistance::Init();
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++){
	    sim_treelets(i,j) = exp(- pow(ed(type_treelet[i],code_treelet[i],
					     type_treelet[j],code_treelet[j]),2)/(options.sigma_inter));
	    if(sim_treelets(i,j) < DBL_EPSILON)
	      sim_treelets(i,j) = 0.0;
	    // cout << sim_treelets(i,j);
	    // if(j!=nb_treelets-1)
	    //   cout <<  ",";
	    // else
	    //   cout << endl;
	  }
	//	sim_treelets.display();
	//	Regul

	CImg<double> eigvals;
	CImg<double> eigvects;
	sim_treelets.symmetric_eigen(eigvals,eigvects);
	cout << "regul ..." << endl;
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++)
	    assert(sim_treelets(i,j) == sim_treelets(j,i));
	if(eigvals[nb_treelets-1] < 0)
	  for (unsigned int i=0; i<nb_treelets; ++i)
	    sim_treelets(i,i) -= eigvals[nb_treelets-1];
	//sim_treelets.symmetric_eigen(eigvals,eigvects);
	//sim_treelets.display();
	((InterTreeletKernel*)kgraph)->setSimTreelets(sim_treelets);
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
      kgraph = new ContractedCycleKernel(options.spectrumKernelContractedCycle,options.sigma);
      if(options.weights_file == NULL){
	if(options.top == 0){
	  vector<string> * vit_list = getSizeContractedCycleTreeletsasVIT(*dataset, options.size_max);
	  ((ContractedCycleKernel*)kgraph)->selectTreelets(vit_list);
	  ((ContractedCycleKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}else if (options.top ==7) {
	  cout << "Selection des treelets d'azote" << endl;
	  //MoleculeGraph util;
	  vector<string> * vit_list = getContractedCycleAzoteAsVIT(*dataset);
	  ((ContractedCycleKernel*)kgraph)->selectTreelets(vit_list);
	  ((ContractedCycleKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
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
    default:
      cerr << "Error : Invalid kernel" << endl;
      exit(1);
    }

  dataset->computeGramMatrix(kgraph, false);
  cout << "Gram Matrix computed" << endl;


  // CV

  // int foldValue=options.foldCV;
  // if(options.foldCV==0)
  //   foldValue=N;

  // int mod=N%foldValue;
  // int tailleFold=N/foldValue;

  // double * cv_bp_exps = new double[N];
  // double * cv_bps = new double[N];
  // double cv_bp_exp_mean = 0.0;
  // double cv_bp_mean = 0.0;

  // int currentMol=0;
  // for(int i=0;i<foldValue;i++)
  //   {
  //     int tailleIter = tailleFold;
  //     if(i<mod)
  // 	tailleIter++;

  //     double* cv_bp_exp =  new double[tailleIter];
  //     Collection** col =  new Collection*[tailleIter];

  //     for(int j=0;j<tailleIter;j++)
  // 	{
  // 	  cv_bp_exp[j] = dataset->getParameter(0);
  // 	  col[j] = dataset->getCollection(0);
  // 	  dataset->delete_first();
  // 	}
      
  //     for(int j=0;j<tailleIter;j++)
  // 	{
  // 	  cv_bps[currentMol] = predict(*dataset,kgraph,options.alpha,col[j],options.regMethode,options.c);
  // 	  cv_bp_exps[currentMol]=cv_bp_exp[j];
  // 	  cv_bp_exp_mean+=cv_bp_exps[currentMol];
  // 	  cv_bp_mean+=cv_bps[currentMol];
  // 	  cout << "CV : Molecule " << currentMol << " : " << cv_bps[currentMol] << " (" << cv_bp_exps[currentMol] << ")" << endl;

  // 	  currentMol++;
  // 	}
      
  //     for(int j=0;j<tailleIter;j++)	  
  // 	dataset->add (col[j], cv_bp_exp[j], kgraph);

  //     delete [] cv_bp_exp;
  //     delete [] col;
  //   }

  // cv_bp_exp_mean /= N;
  // cv_bp_mean /= N;


  // CV (no fold)
  double * cv_bp_exps = new double[N];
  double * cv_bps = new double[N];
  double cv_bp_exp_mean = 0.0;
  double cv_bp_mean = 0.0;

  // for(int i=0;i<N;i++)
  //   {
  //     cv_bps[i] = predict(*dataset,kgraph,options.alpha,dataset->getCollection(i),options.regMethode,options.c);
  //     cv_bp_exps[i]=dataset->getParameter(i);
  //     cv_bp_exp_mean+=cv_bp_exps[i];
  //     cv_bp_mean+=cv_bps[i];
  //     cout << "CV : Molecule " << i << " : " << cv_bps[i] << " (" << cv_bp_exps[i] << ")" << endl;
  //   }

  cv_bp_exp_mean /= N;
  cv_bp_mean /= N;


  //Test molecules prediction
  list<int>::iterator it = to_test.begin();
  vector<Collection*>::iterator it_col = to_test_Collections.begin();    
  list<double>::iterator it_param = to_test_Parameters.begin();
    
  for(;it != to_test.end();++it, ++it_col, ++it_param){
    clock_t predict_start,predict_end;
    struct timeval predict_s_start, predict_s_end;
    predict_start = gettimeofday(&predict_s_start,NULL);
    double bp = predict(*dataset, kgraph, options.alpha,*it_col, options.regMethode, options.c, options.svm_epsilon);
    predict_end = gettimeofday(&predict_s_end,NULL);//times(&s_end);
    mean_prediction_time += (1000000*predict_s_end.tv_sec + predict_s_end.tv_usec) - 
      (1000000*predict_s_start.tv_sec + predict_s_start.tv_usec);
    bps[nb_molecules]=bp;
    double bp_exp = *it_param;
    double err =  bp_exp - bp;
    err_quad +=  err * err;
    err_moy += abs(bp_exp - bp);
    bp_mean += bp;
    cout << "Molecule " << nb_molecules << " : " << bp << " (" << bp_exp << ")" << endl;
    nb_molecules ++;
  }
 
  cout<< "Time predict "<<mean_prediction_time<<endl;


  bp_mean /= N2;
  std_dev = sqrt(err_quad / N2); 
  err_moy /= N2;

  double q_den = 0.0;
  for(int i =0;i < N; i++){
    q_den += pow(cv_bp_exps[i] - cv_bp_exp_mean,2);
  }

  double q_num =0.0;
  for(int i =0;i < N; i++){
    q_num += pow(cv_bp_exps[i] - cv_bps[i],2);
  }

  double q=1.0 - (q_num / q_den);

  double r_den = 0.0;
  for(int i =0;i < N2; i++){
    r_den += pow(bp_exps[i] - bp_exp_mean,2);
  }

  double r_den2 = 0.0;
  for(int i =0;i < N2; i++){
    r_den2 += pow(bps[i] - bp_mean,2);
  }
  double r_num =0.0;
  for(int i =0;i < N2; i++){
    r_num += (bps[i] - bp_mean)*(bp_exps[i] - bp_exp_mean);
  }

  double r2 = r_num / (sqrt(r_den*r_den2));

  //Print Results
  cout << endl;
  cout << "Average error : " << err_moy  << endl;
  cout << "Standard Deviation : " << std_dev << endl;
  cout<<endl;
  cout << "q² : " << q << endl;
  cout << "R : " << r2 << endl;
  if(q>0 && r2>0)
    {
      cout << "Norme : " << sqrt(q*q+r2*r2) << endl;
    }

  delete edit;
  return 0;
}

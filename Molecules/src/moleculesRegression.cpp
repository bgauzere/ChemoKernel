/*
 * @file moleculesRegression.cpp
 *
 * This program performs a Kernel Ridge regression on molecules datasets.
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-10)
 */


#include <iostream>
#include <fstream>

#include <algorithm>
#include <cfloat>
#include "MoleculeGraph.h"
#include "GraphKernel.h"
#include "TrailKernel.h"
#include "ShapeBagTrails.h"
#include "TreeletType.h"
#include <pandore.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sstream>
#include <assert.h>
#include "KMean.h"
#include "KWMean.h"
#include "RandomWalkKernel.h"
#include "KRegression.h"
#include "GraphletCountKernel.h"
#include "TreeletKernel.h"
#include "CombinedKernel.h"
#include "RandomWalkKernel.h"
#include "Kashima.h"
#include "TanimotoKernel.h"
#include "MinMaxKernel.h"
#include "KernelRidge.h"
#include "KernelRidgeTopN.h"
#include "WeightedAverage.hpp"
#include "SVM.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "LaplacianKernelOriginal.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistanceV2.h"
#include "KEditDistance.h"
#include "KPrototypeEditDistance.h"
#include <libgen.h>
#include "utils.h"
#include "string_utils.h"
#include "TreeletCycleKernel.h"
#include "InterTreeletKernel.h"
#include "GaussianKernel.h"
#include "InnerProductKernel.h"
#include "TreeletEditDistance.h"
#include "TreeletIntersectionKernel.h"
#include "TreeletGaussianKernel.h"
#include "TreeletInnerProductKernel.h"
#include "TreeletBinaryKernel.h"
#include "TreeletRandomKernel.h"
#include "TreeletCompleteGaussianKernel.h"
extern "C" {
#include "gnuplot_i.h"
}

using namespace std;
using namespace pandore;

#define LIMIT_LOW_FREQUENCY 140
bool always_true(double i)
{
  return true;
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


vector<string> * getAcyclicVIT()
{
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];  
  
  vit_list[0].push_back("C");
  vit_list[0].push_back("S");

  vit_list[1].push_back("C1O");
  vit_list[1].push_back("O1O");
  vit_list[1].push_back("C1S");

  vit_list[2].push_back("O1C1O");
  vit_list[2].push_back("C1O1O");

  vit_list[3].push_back("C1C1C1C");
  vit_list[3].push_back("C1C1C1S");
  vit_list[3].push_back("C1S1S1C");
  vit_list[3].push_back("S1C1C1S");

  vit_list[4].push_back("O1C1C1C1O");
  vit_list[4].push_back("C1C1C1C1O");
  vit_list[4].push_back("C1C1S1C1S");
  vit_list[4].push_back("C1C1O1C1C");
  vit_list[4].push_back("C1S1C1C1S");
  vit_list[4].push_back("C1C1C1O1O");
  
  vit_list[5].push_back("C1O1C1C1C1O");
  vit_list[5].push_back("C1C1C1C1S1C");
  
  vit_list[6].push_back("C1C1C1S");
  vit_list[6].push_back("C1C1C1O");
  vit_list[6].push_back("C1C1C1C");
  vit_list[6].push_back("C1C1S1S");

  vit_list[9].push_back("C1C1C1C1C1C");  
  //  vit_list[9].push_back("C1C1C1C1C1S");  

  vit_list[10].push_back("O1C C1C1C1C1C");

  return vit_list;
}

vector<string> * getVilleminVIT()
{
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];  
  
  vit_list[0].push_back("C");
  vit_list[0].push_back("O");  
  vit_list[0].push_back("S");
  
  vit_list[1].push_back("C1C");
  vit_list[1].push_back("C1S");
  vit_list[1].push_back("C1O");
  vit_list[1].push_back("S1S");
  vit_list[1].push_back("O1O");
  
  vit_list[2].push_back("C1C1S");
  vit_list[2].push_back("C1C1C");
  vit_list[2].push_back("C1S1C");
  vit_list[2].push_back("C1O1C");
  vit_list[2].push_back("S1C1S");
  vit_list[2].push_back("C1C1O");
  vit_list[2].push_back("O1C1O");

  vit_list[3].push_back("C1C1C1C");
  
  vit_list[4].push_back("C1C1C1C1C");

  vit_list[6].push_back("C1C1C1C");

  vit_list[7].push_back("C1S1C1C1C");

  vit_list[8].push_back("C1C1C1C1C");
  
  return vit_list;
}


bool comp (pair<string, double>  l, pair<string, double>  r)
{ 
  return (fabs(l.second) > fabs(r.second)); 
}

vector<string> * getVIT(MoleculesDataset dataset,int top)
{
  dataset.computeLabeledGraphletCorrelation();
  map<string, double,bool (*)(string, string)> * correlation = dataset.getLabeledCorrelation();
  vector< pair<string, double> > sort_correlations;
  
  for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
    {
      stringstream out;
      out << "-G" << j;
      string code_treelet= out.str();
      map<string, double,bool (*)(string, string)>::iterator it = correlation[j].begin();
      for(;it != correlation[j].end();it ++)
	{
	  string nouveau_code(it->first);
	  nouveau_code.append(code_treelet);
	  sort_correlations.push_back(pair<string, double>(nouveau_code, it->second));
	}
    }
  //sort_correlations classé par corrélation décroissante
  sort(sort_correlations.begin(), sort_correlations.end(), comp);
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  for(int i=0;i<top;i++)
    {
      cout << sort_correlations[i].first << endl;
      string tmp;
      int start = sort_correlations[i].first.length()-2;
      tmp = sort_correlations[i].first.substr(start,2);
      stringstream parse;
      parse << tmp;
      int treelet_type;
      parse >> treelet_type;
      int length_code = sort_correlations[i].first.length() - 4;
      string code = sort_correlations[i].first.substr(0,length_code);
      vit_list[treelet_type].push_back(code);
    }
  
  return vit_list;
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
    //  while(nb_vit != 56) // >0
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
		  ;
		}
	      delete [] tmp_vit_list;
	    }
	  
	}
      // for(int treelet_type=0;treelet_type<SIZE_SPECTRUM;treelet_type++)
      // 	for(int i=0;i<candidate_vit_list[treelet_type].size();i++)
      // 	  cout << treelet_type << candidate_vit_list[treelet_type][i] << endl;
      // cout << endl << nb_vit << " " << std_dev_current << endl;
      //Ici, on a un nouvel ensemble de vit + un nouveau ecart type 
      nb_vit --;
    }
  return  vit_list;
}


vector<string> * getStepWiseVIT (MoleculesDataset & dataset,  TreeletKernel* kgraph, double lambda)
{	
  //XXX: On calcule la corrélation en considérant ** TOUTE ** la base, ce qui est faux
  //     puisqu'on est pas sensé connaitre la temp réelle de la molécule à tester.
  //dataset.computeLabeledGraphletCorrelation();

  dataset.computeLabeledSpectrums();
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  int N = dataset.size();
  double std_dev = DBL_MAX;
  double std_dev_current = 80;
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  // for(int i=0;i<SIZE_SPECTRUM;i++)
  //   vit_list[i] = vector<char *>;
  vector<string> * candidate_vit_list = new vector<string>[SIZE_SPECTRUM];
  //XXX: Vieux Hack de sioux pour eviter le produit scalaire de deux trucs vides 
  candidate_vit_list[0].push_back("C");  

  int nb_regressions = 0;
  int nb_vit = 0;
  while(std_dev_current < std_dev)
    { //On s'arrete quand on progresse plus, (évolution a plotter ds un 2nd temps)
      nb_vit ++;
      std_dev = std_dev_current;
      std_dev_current = DBL_MAX;
      for(int i=0;i<SIZE_SPECTRUM;i++)
	vit_list[i] = vector<string>(candidate_vit_list[i]);
      
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


vector<string> * getNStepWiseVIT (MoleculesDataset & dataset,  TreeletKernel* kgraph, 
				  double lambda, int frequent_threshold, int step)
{	
  //Forward Selection 3par3 
  dataset.computeLabeledSpectrums();
  //On supprime les non fréquents
  int nb_frequent_treelets = 0;
  treelet_spectrum * distribution = dataset.getTreeletDistribution(&always_true);
  for(int i =0; i < SIZE_SPECTRUM; i++)
    {
      vector <string> to_erase;
      treelet_spectrum::iterator it = distribution[i].begin();
      for(;it != distribution[i].end();it ++)
	if(it->second < frequent_threshold)
	  to_erase.push_back(it->first);
	else
	  nb_frequent_treelets ++;
      for(unsigned int j=0;j<to_erase.size();j++)
	distribution[i].erase(distribution[i].find(to_erase[j]));
    }
  distribution[0].erase("C");  
  nb_frequent_treelets --;
    
  
  //Verif
  // for(int i =0; i < SIZE_SPECTRUM; i++)
  //   {
  //     MoleculeGraph::spectrum_map::iterator it = distribution[i].begin();
  //     for(;it != distribution[i].end();it ++)
  // 	cout << it->second << endl;
  //   }
  int N = dataset.size();
  double std_dev = DBL_MAX;
  double std_dev_current = 80;
  //Liste élue
  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
  //Liste courante de  treelets a tester
  vector<string> * candidate_vit_list = new vector<string>[SIZE_SPECTRUM];
  //XXX: Vieux Hack de sioux pour eviter le produit scalaire de deux trucs vides 
  candidate_vit_list[0].push_back("C");  
  
  int nb_regressions = 0;
  int max_vit = nb_frequent_treelets; //-1 : "C"
  cout << nb_frequent_treelets << endl;
  int nb_vit = 1;
  //while(std_dev_current < std_dev) //On s'arrete quand on progresse plus
  while(nb_vit < max_vit) //On s'arrete quand y a plus rien
    { 
      nb_vit += step;
      std_dev = std_dev_current;
      std_dev_current = DBL_MAX;
      //Maj des gagnants
      for(int i=0;i<SIZE_SPECTRUM;i++)
	{
	  vit_list[i] = vector<string>(candidate_vit_list[i]);
	  //Affichage des gagnants
	  for(unsigned int j = 0;j<vit_list[i].size();j++)
	    {
	      cout << "G" << i << "-" << vit_list[i][j] << endl;
	      //Sortie de distribution
	      distribution[i].erase(vit_list[i][j]);
	    }
	}
      
      //Construction des permutations possibles
      string* treelets_to_test = new string[nb_frequent_treelets];
      int * map_code_struture = new int [nb_frequent_treelets];/*permet de retrouver la 
								 structure associée au code de l'indice*/
      int j = 0;
      for(int i =0; i < SIZE_SPECTRUM; i++)
	{
	  treelet_spectrum::iterator it = distribution[i].begin();
	  for(;it != distribution[i].end();it ++)
	    {
	      treelets_to_test[j] = it->first;
	      map_code_struture[j] = i;
	      j++;
	    }
	}
        string ** permutations_treelets;
	int ** permutations_positions;
	cout << "Calcul des permuts"<< endl;
	string null_string("");
	int nb_candidats = utils::get_perm(treelets_to_test, j, 
					   permutations_treelets, permutations_positions, step,  null_string); 
	/***********************/
	cout << nb_candidats << endl;
	/***********************/
	
	for(int j=0;j<nb_candidats;++j)//Parcours de toutes les permutation possibles
	  {
	    //Regression avec permutations
	    //Creation de la vit_list candidate
	    vector<string> * tmp_vit_list = new vector<string>[SIZE_SPECTRUM];
	    for(int i=0;i < SIZE_SPECTRUM;i++)
	      tmp_vit_list[i] = vector<string>(vit_list[i]);
	    
	    //Rajout des treelets candidats
	    for(int n=0;n<step;n++)
	      {
		tmp_vit_list[map_code_struture[permutations_positions[j][n]]].push_back(permutations_treelets[j][n]);
	      }
	    kgraph->selectTreelets(tmp_vit_list); //XXX:Rajouter les poids
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
	    cout << "Régression " << nb_regressions << "/" << nb_candidats  << endl;
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
	
	cout << nb_vit << " " << std_dev_current << endl;
	//Ici, on a un nouveau vit + un nouveau ecart type 
    }
  return  vit_list;
}


void leaveOneOut (MoleculesDataset & dataset, GraphKernel* kgraph, double lambda,bool quietMode, int nb_neighbours,int methode,double c,double epsilon)
{	
  KRegression* kr;
  switch(methode)
    {
    case 1: kr=new KernelRidge(kgraph, &dataset, lambda); break;
    case 2: kr=new KernelRidgeTopN(kgraph, &dataset, lambda, nb_neighbours); break;
    case 3: kr= new WeightedAverage(kgraph, &dataset); break;
    default: kr=new KernelRidge(kgraph, &dataset, lambda); break;
    }

  double err_moy = 0; 
	
  int N = dataset.size();

  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_exp_mean = 0.0;
  double bp_mean = 0.0;
  double * bp_exps = new double[N];
  double * bps = new double[N];

  for (int i=0; i<N; ++i)
    {
      double bp_exp = dataset.getParameter(0);
      Collection* col = dataset.getCollection(0);
      dataset.delete_first();

      double bp;
      if(methode==4)
	{
	  SVM svm (dataset, kgraph, SVM::EPSILON_SVR,c,0.5,epsilon);
	  bp=svm.predict(col);
	}
      else
	bp = (*kr) (col);
      
      //trop dependant de la temp d ebullition (grosse erreur sur temp proche de 0)
      // cout << "Molecule " << i << " : " << bp << " (" << bp_exp << ")";
      // cout << " --> " << abs(bp_exp - bp)/bp_exp*100 << "% (lambda = " << lambda << ")" << endl;
      // err_moy += (bp_exp == 0)?bp:abs(abs(bp_exp - bp)/bp_exp*100);
		
      double err = bp_exp - bp;
      err_quad += err * err;
		
      if(!quietMode){
	cout << "Molecule " << i << " : " << bp << " (" << bp_exp << ")";
	cout << "\t --> \t" << abs(err)  << "° (lambda = " << lambda << ")" << 
	  dataset.getMoleculeFilename(col->GETVALUE("id",Long))<< endl;
      }
      err_moy += abs(bp_exp - bp);
      bp_exp_mean += bp_exp; 
      bp_mean += bp; 
      bp_exps[i] = bp_exp;	
      bps[i] = bp;
      dataset.add (col, bp_exp, kgraph);
    }

  bp_mean /= N;
  bp_exp_mean /= N;
  std_dev = sqrt(err_quad / N); 
  err_moy /= N;

  double r_den = 0.0;
  for(int i =0;i < N; i++){
    r_den += pow(bp_exps[i] - bp_exp_mean,2);
  }
  double r = sqrt(1- (err_quad) / (r_den));

  cout << endl;
  cout << "Average error : " << err_moy  << endl;
  cout << "Standard Deviation : " << std_dev << endl;
  cout << "R : " << r << endl;
  cout<<endl;


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
  double mu;
  bool normalize;
  KernelType spectrumKernel;
  KernelType spectrumKernelCycle;
  int top;
  int nb_neighbours;
  int step;
  char * weights_file;
  KPrototypeEditDistance::PrototypeSelectionType prototypeSelection;
  int regMethode;
  int phiPow;
  double c;
  int labelCycles;
  bool noUnique;
  double delta;
  int nbAtomsMax;
  double epsilon;
  int epsiPropSD;
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
  cerr << "\t5 -- Random Walks Kernel" << endl;
  cerr << "\t6 -- Treelet Kernel" << endl;
  cerr << "\t7 -- Tanimoto Kernel on Treelets" << endl;
  cerr << "\t8 -- MinMax Kernel on Treelets" << endl;
  cerr << "\t9 -- Test Graph Laplacian Kernel" << endl;
  cerr << "\t10 -- Edit Distance Prototype Selection" << endl;
  cerr << "\t11 -- Treelet Cycle Kernel" << endl;
  cerr << "\t13 -- Inter Treelet Kernel - Approximate Edit Distance" << endl;
  cerr << "\t14 -- Inter Treelet Kernel - Laplacian Kernel" << endl;
  cerr << "\t15 -- Inter Treelet Kernel - Exact Connexe Edit Distance (Regul)" << endl;
  

  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Binary Kernel" << endl;
  cerr << "\t4 -- Random Kernel" << endl;
  cerr << "\t5 -- Complete Gaussian Kernel" << endl;
  cerr << "\t1 -- V1 for acyclic - Parcours par couche pour cas prenant en compte les cycles" << endl;
  cerr << "\t2 -- V2 for acyclic - Parcours par branche pour cas prenant en compte les cycles" << endl;

  cerr << "-P method : Set the method to select the prototypes" << endl;
  cerr << "-r nb proto : using -k 10, set the number of prototypes" << endl;

  cerr << "-a alpha : Set the parameter of the regression (default 1.0)" << endl;
  cerr << "-s sigma : Set the sigma parameter of the Gaussian Kernel & Graph Laplacian kernel & Graphlet Count Kernel (default 2.0)" << endl;
  cerr << "-r regularization : Set the regularization type of the Graph Laplacian Kernel (default 0)" << endl;
  cerr << "-l lambda : Set the lambda regularization parameter of the Graph Laplacian kernel (default 1.0)" << endl;
  cerr << "-f keptTrails : Set the percent of kept trails in bag of trails kernels (KMean, KWMean) (default 1.0)" << endl;
  cerr << "-m trailsLength : Set the maximum length of trails in bag of trails kernels (KMean, KWMean) (default 15)" << endl;
  cerr << "-e editdistance : Set the Edit Distance Costs: " << endl;
  cerr << "\t0 -- V1" << endl;
  cerr << "\t1 -- V2" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-N : Normalize Treelet Kernel" << endl;
  cerr << "-t n : Sélectionne les n treelets les plus représentatifs " << endl;
  cerr << "-b n : Sélectionne les n molecules les plus proches pour la régression" << endl;
  cerr << "-x n : Fait une stepwise selection n par n" << endl;
  cerr << "-w weigths : Utilise les poids weights, combinée avec -t -7" << endl;
  cerr << "-reg regression : Choisis différentes méthodes de régression" << endl;
  cerr << "\t1 -- Kernel Ridge" << endl;
  cerr << "\t2 -- Kernel Ridge Top N" << endl;
  cerr << "\t3 -- Weighted Average of Known Values" << endl;
  cerr << "\t4 -- SVM" << endl;
  cerr << "-c C SVM : Paramètre C du SVM" << endl; 
  cerr << "-eps epsilon : Paramètre Epsilon du SVM " << endl;
  cerr << "-epsSD : Epsilon est prit comme proportionnel (à l'inverse de la valeur passer par -eps ) à l'écart type du jeu de données" << endl;
  cerr << "-lC Modifications of label"<< endl;
  cerr << "\t0 -- No modifications" << endl;
  cerr << "\t1 -- Number of relevant cycles add to the label of an atom" << endl;
  cerr << "\t2 -- Size of each relevant cycles add to the label of an atom" << endl;
  cerr << "-u supprimer les treelets présents dans une seul molécule" << endl;
  cerr << "\t -nb nombre maximum d'atome par treelets" << endl;


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
  options->normalize = false;
  options->spectrumKernel = InnerProductKernelType;
  options->spectrumKernelCycle = InnerProductKernelType;
  options->top = 0;
  options->nb_neighbours = 10;
  options->step = 3;
  options->weights_file = NULL;
  options->regMethode=1;
  options->phiPow = 2;
  options->c = 1.0;
  options->labelCycles=0;
  options->noUnique=false;
  options->delta = 0.5;
  options->nbAtomsMax=SIZE_MAX;
  options->epsilon=1;
  options->epsiPropSD=0;


  int i=3;
  while (i<argc)
    {
     if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->delta = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-epsSD", 6) == 0)
	{
	  options->epsiPropSD=1;
	  i+=1;
	}
      else if (strncmp(argv[i], "-eps", 4) == 0)
	{
	  options->epsilon = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->kernel = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-lC", 3) == 0)
	{
	  options->labelCycles = atoi(argv[i+1]);
	  i+=2;
	} 
      else if (strncmp(argv[i], "-reg", 4) == 0)
	{
	  options->regMethode = atoi(argv[i+1]);
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
      else if (strncmp(argv[i], "-nb", 3) == 0)
	{
	  options->nbAtomsMax = atoi(argv[i+1]);
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
      else if (strncmp(argv[i], "-p", 2) == 0)
	{
	  options->phiPow = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-u", 2) == 0)
	{
	  options->noUnique = true;
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
	
  ShapeBagTrails* sbt = new ShapeBagTrails (options.trailsLength, options.keptTrails, 0.0f, false, false, 0);
	
  GraphKernel* kgraph = NULL;
  GraphEditDistance* edit = NULL;
  if(options.editdistance == 0){
    edit = new MoleculeGraphEditDistance;
  }else{
    //V2
    edit = new MoleculeGraphEditDistanceV2;
  }
	
  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.trainset_path, options.dataset_file);
  //dataset.setAbsParameter();

  if(options.epsiPropSD)
    options.epsilon=dataset.standardDeviation()/options.epsilon;

  if(options.regMethode==4)
    cout<<"Epsilon : "<<options.epsilon<<endl;

  double * coeffs;
  long diff_2,diff;	
  clock_t start,end;
  struct timeval s_start, s_end;
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
      
      start = gettimeofday(&s_start,NULL);
      kgraph = new LaplacianKernel (edit, dataset, options.sigma, options.regularization, options.lambda);
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff_2 = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Laplacian Kernel Initialisation : " << diff_2 << endl;
      cout << "Poids :" << endl;
      //((LaplacianKernel *)kgraph)->printWeightMatrix();
      break;
    case 3:
      kgraph = new KEditDistance (edit, options.sigma);
      break;	
    case 4:  
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
      kgraph = new RandomWalkKernel(new Kashima(),options.lambda,1.0);
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
	case CompleteGaussianKernelType:
	  cout << "Complete Gaussian Kernel";
	  break;
	case JaccardKernelType:
	  cout << "Jaccard Kernel";
	  break;

	default:
	  //devrait pas arriver
	  break;
	}
      cout << endl;
      
      start = gettimeofday(&s_start,NULL);
      dataset.computeLabeledSpectrums();
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff_2 = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Spectrum Computation Time : " << diff_2 << endl;
      dataset.getTreeletDistribution(&always_true);
      cout << "Nb Treelets = " << dataset.getNbTreelets() << endl;
      if(options.normalize)
	dataset.normalizeLabeledSpectrums();
      kgraph = new TreeletKernel(options.spectrumKernel,options.sigma);
      if(options.top > 0)
	{
	  vector<string> * vit_list = getVIT(dataset,options.top);
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      if(options.top == 0)
	{
	  vector<string> * vit_list = getAllAsVIT(dataset);
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}

      else if (options.top == -1)
	{
	  vector<string> * vit_list = getVilleminVIT();
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -2)
	{
	  vector<string> * vit_list = getStepWiseVIT(dataset,(TreeletKernel*)kgraph,options.alpha);
	  cout << "VIT List : " << endl;
	  for(int i=0;i<SIZE_SPECTRUM;i++)
	    for(unsigned int j=0;j<vit_list[i].size();j++)
	      cout << "G" << i << "-" << vit_list[i][j] << endl;
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -3)
	{
	  vector<string> * vit_list = getAcyclicVIT();
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -4)
	{
	  vector<string> * vit_list = getAcyclicVIT_2();
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  map<string, double,bool (*)(string, string)>* weights = vit_to_weight(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(weights);
	  
	}
      else if (options.top == -5)
	{
	  vector<string> * vit_list = getNStepWiseVIT(dataset,
						      (TreeletKernel*)kgraph,
						      options.alpha, 
						      LIMIT_LOW_FREQUENCY,options.step);
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -6)
	{
	  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
	  vit_list[0].push_back("C");
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -7)
	{
	  vector<string> * vit_list = new vector<string>[SIZE_SPECTRUM];
	  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
	  map<string, double,bool (*)(string, string)> * weights = 
	    new map<string, double,bool (*)(string, string)>[SIZE_SPECTRUM];
	  for(int k=0;k<SIZE_SPECTRUM;k++)
	    {
	      weights[k] = map<string, double, bool (*)(string, string)> (fn_pt); 
	    }
  
	  //lecture du fichier
	  ifstream f (options.weights_file);
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
	      parse_line << line.substr(pos_sep_code_weight+1,line.npos-pos_sep_code_weight-1);
	      parse_line >> weight;
	      weights[treelet_type].insert(pair<string, double>(code,weight));
	      vit_list[treelet_type].push_back(code);
	    }
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(weights);
	} 
      else if (options.top == -8)//Backward Elimination
	{
	  vector<string> * vit_list = getBackwardStepWiseVIT(dataset,(TreeletKernel*)kgraph,options.alpha);
	  cout << "VIT List : " << endl;
	  for(int i=0;i<SIZE_SPECTRUM;i++)
	    for(unsigned int j=0;j<vit_list[i].size();j++)
	      cout << "G" << i << "-" << vit_list[i][j] << endl;
	  ((TreeletKernel*)kgraph)->selectTreelets(vit_list);
	  ((TreeletKernel*)kgraph)->weightTreelets(vit_to_weight(vit_list));
	}
      else if (options.top == -9)//Backward Elimination
	{
	  MoleculeGraph util;
	  util.initTable ();
	  vector<string> * vit_list = getBackwardStepWiseVIT(dataset,(TreeletKernel*)kgraph,options.alpha);
	  cout << "VIT List for All Molecules : " << endl;
	  for(int i=0;i<SIZE_SPECTRUM;i++)
	    for(unsigned int j=0;j<vit_list[i].size();j++)
	      printf("G%.2d-%s\n",i,util.translateTreeletCode(string(vit_list[i][j])).c_str());
		     
	  int N = dataset.size();
	  for (int i=0; i<N; ++i)
	    {
	      cout << "VIT List without Molecule " << i << " : "  << endl;
	      double bp_exp = dataset.getParameter(0);
	      Collection* col = dataset.getCollection(0);
	      dataset.delete_first();
	      vit_list = getBackwardStepWiseVIT(dataset,(TreeletKernel*)kgraph,options.alpha);
	      for(int i=0;i<SIZE_SPECTRUM;i++)
		for(unsigned int j=0;j<vit_list[i].size();j++)
		  printf("G%.2d-%s\n",i,util.translateTreeletCode(string(vit_list[i][j])).c_str());
	      dataset.add (col, bp_exp, kgraph);
	    }
	  exit(EXIT_SUCCESS);
	}
      
      break;
	    
    case 7:  
      cout << "Kernel : Tanimoto Kernel on Treelets" << endl;
      dataset.computeLabeledSpectrums();
      kgraph = new TanimotoKernel();
      break;
    case 8:  
      cout << "Kernel : MinMax Kernel on Treelets" << endl;
      dataset.computeLabeledSpectrums();
      kgraph = new MinMaxKernel();
      break;
    case 9:
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
      
      kgraph = new LaplacianKernelOriginal (edit, dataset, options.sigma, options.regularization, options.lambda);
      break;
    case 10:
      cout << "Kernel : Prototype Edit Distance" << endl;
      start = gettimeofday(&s_start,NULL);
      kgraph = new KPrototypeEditDistance(dataset, options.regularization,options.prototypeSelection);
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff_2 = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Prototype Selection Embedeed Time : " << diff_2 << endl;
 
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
      dataset.computeRelevantCycles();
      start = gettimeofday(&s_start,NULL);
      dataset.computeCyclesSpectrums();
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Cycle Spectrum Computation Time : " << diff << endl;
      
      kgraph = new CycleKernel(options.spectrumKernelCycle,options.sigma);
      
      if(options.weights_file == NULL){
	((CycleKernel*)kgraph)->selectTreelets(getCyclesAllAsVIT(dataset));
	((CycleKernel*)kgraph)->weightTreelets(vit_to_weight(getCyclesAllAsVIT(dataset)));
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
      dataset.getCycleTreeletDistribution(&always_true);
      cout << "Nb Cycle Treelets : " << dataset.getNbCyclesTreelets() << endl;
      break;
    case 13:      
      {
	dataset.computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
	int nb_treelets =  dataset.getNbTreelets();
	cout << nb_treelets<< endl;
	dataset.computeSpecialVectors(&distrib, nb_treelets);
	cout << "Special Vectors computed " << endl;
	
	GraphKernel * kinter_treelets = new KEditDistance (edit, 1.0);
	Kernel * kinter_graphs = new GaussianKernel(1.0);
	kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
					&distrib,
					nb_treelets);
      }
      break;
    case 14:      
      {
	dataset.computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
	int nb_treelets =  dataset.getNbTreelets();
	cout << nb_treelets<< endl;
	dataset.computeSpecialVectors(&distrib, nb_treelets);
	cout << "Special Vectors computed " << endl;
	
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
	// GraphKernel * kinter_treelets = new KEditDistance (edit, 1.0);
	GraphKernel * kinter_treelets = new LaplacianKernelOriginal (edit, dataset_treelets,
								     options.sigma,
								     options.regularization,
								     options.lambda);
	cout << "Laplacian Computed" << endl;
	Kernel * kinter_graphs = new GaussianKernel(2.0);
	cout << "TreeletKernel Initialized" << endl;
	kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
					col_treelets,
					nb_treelets);
      }
      break;
    case 15:      
      {
	dataset.computeLabeledSpectrums();
	treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
	int nb_treelets =  dataset.getNbTreelets();
	cout << nb_treelets<< endl;
	dataset.computeSpecialVectors(&distrib, nb_treelets);
	cout << "Special Vectors computed " << endl;
	
	Kernel * kinter_graphs = new GaussianKernel(options.sigma);
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
	
	TreeletEditDistance ed(1,7);	
	TreeletEditDistance::Init();
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++){
	    sim_treelets(i,j) = exp(- pow(ed(type_treelet[i],code_treelet[i],
					     type_treelet[j],code_treelet[j]),2)/(5));
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
	for (unsigned int i=0; i<nb_treelets; ++i)
	  cout << eigvals[i] << endl;
	cout << "regul ..." << endl;
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++)
	    assert(sim_treelets(i,j) == sim_treelets(j,i));
	if(eigvals[nb_treelets-1] < 0)
	  for (unsigned int i=0; i<nb_treelets; ++i)
	    sim_treelets(i,i) -= eigvals[nb_treelets-1];
	sim_treelets.symmetric_eigen(eigvals,eigvects);
	for (unsigned int i=0; i<nb_treelets; ++i)
	  cout << eigvals[i] << endl;
	//sim_treelets.display();
	((InterTreeletKernel*)kgraph)->setSimTreelets(sim_treelets);
      }
      break;
    default:
      cerr << "Error : Invalid kernel" << endl;
      exit(1);
    }
  if(!options.quietMode){
    cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
    cout << "Number of molecules : " << dataset.size() << endl;
    cout << "Alpha : " << options.alpha << endl;
    cout << endl;
  }
	
  // We compute the gram matrix of the training set
  cout << "Gram Matrix Computation ...";
  cout.flush();
  start = gettimeofday(&s_start,NULL);
  dataset.computeGramMatrix(kgraph, false);
  end = gettimeofday(&s_end,NULL);//times(&s_end);
  cout << "Done" << endl;
  diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
    (1000000*s_start.tv_sec + s_start.tv_usec);
  cout.precision(DBL_DIG);
  cout << "Gram Matrix Computation Time : " << diff << endl;
  bool isPD = dataset.isGramMatrixPD();
  cout << "Gram Matrix Positive definite ? ! " << isPD << endl;
  // if(!isPD)
  //   dataset.regularizeGramMatrix();
  
  dataset.showGramMatrixMatlab ("gram.mat");

  
  //dataset.showIdenticalLine();

  //cout << "Gram Matrix Positive definite ? ! " << dataset.isGramMatrixPD() << endl;
  cimg_library::CImg<double> gram = dataset.getGramMatrix (false);
  
  // int N = dataset.size();
  // ofstream outfile_gram ("gram.mat");
  // outfile_gram.precision(9);
  // //outfile_gram << N << " " << N << endl;
  // outfile_gram.flush();
  // for (int i=0; i<N; ++i)
  //   {
  //     for (int j=0; j<N; ++j)
  // 	{
	  
  // 	  outfile_gram << gram(i,j);
  // 	  outfile_gram.flush();
  // 	  if(j!=N-1)
  // 	    outfile_gram.flush() << ",";
  // 	}
  //     outfile_gram << endl;
  //   }
  // outfile_gram.close();
  // cout << "gram.mat Ok."<< endl;

  // int mol = 162;
  // for(int i=0;i<dataset.size();i++){
  //   cout << gram(mol,i) << " = (" << basename(dataset.getMoleculeFilename(mol)) << "(" << mol<< ")" << "." 
  // 	 << basename(dataset.getMoleculeFilename(i)) << "(" << i<< ")" << ")"<< endl;
  // }
  // for(unsigned int i=0;i<dataset.size();i++)
  //   cout << gram(i,i) << endl;
  // for(unsigned int i=0;i<dataset.size();i++)
  //   for(unsigned int j=0;j<dataset.size();j++)
  //     if (gram(i,j) == 1 && (i!=j))
  // 	cout << basename(dataset.getMoleculeFilename(i)) << "(" << i << ")" << " = " 
  // 	     << basename(dataset.getMoleculeFilename(j)) << "(" << j << ")" << endl;
 
  //dataset.normalizeParams(0);
  //cout << "Gram SVM Format " << endl;
  // dataset.showGramMatrix();
  cout << "Gram Raw Format " << endl;
   dataset.showGramMatrixRaw();
	
  // Start the leave-one-out
  //cout << options.nb_neighbours << endl;
   leaveOneOut(dataset, kgraph, options.alpha,options.quietMode, options.nb_neighbours,options.regMethode,options.c,options.epsilon);
  
  // Destruction of the training set
  dataset.destroy();
	
  delete edit;
  delete kgraph;
  return 0;
}


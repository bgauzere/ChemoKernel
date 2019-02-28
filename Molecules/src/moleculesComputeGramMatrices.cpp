
#include <iostream>
#include <sstream>
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
#include "CImg.h"
#include <fstream>
#include <cassert>

using namespace std;
using namespace pandore;
using namespace cimg_library;

void showUsage ()
{
  cerr << "Usage : moleculesClassification trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-o id : File Identifier replaced in gram_[id].mat, treelets_[id].mat, y_[id].mat and treelets_cool_[id].mat" << endl;
  cerr << "-k kernel:" << endl;
  cerr << "\t0 -- Treelet Kernel" << endl;
  cerr << "\t1 -- Cycle Kernel" << endl;
  cerr << "\t2 -- StereoKernel"<<endl;

  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Polynomial Kernel" << endl;
  cerr << "\t4 -- Complete Gaussian Kernel" << endl;

  
  cerr << "-s sigma : Set the parameter for gaussian kernel or power for polynomial kernel" << endl;
  cerr << "-N : Normalize Kernel" << endl;

  cerr << "-gs : One gram matrix by size of extensions" << endl;

}

//XXX: Code a factoriser
bool always_true(double unused)
{
  return true;
}

double polynomialKernel(double t1,double t2, double power){
  return pow((t1*t2+1),power);
}
double innerKernel(double t1,double t2, double unused){
  return t1*t2;
}
double gaussianKernel(double t1,double t2, double sigma){
  return exp(-((t1-t2)*(t1-t2))/(sigma*sigma));
}
double intersectionKernel(double t1,double t2, double unused){
  return min(t1,t2);
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
//Fin code a factoriser

struct Options
{
  char* trainset_path;
  char* dataset_file;
  int kernel;
  int kernel_type;
  double sigma;
  bool normalize;
  bool quietMode;
  char * output_id;
  bool oneGramBySize;
  int nbAtomsMax;
};

void readOptions (int argc, char** argv, Options* options)
{
  if (argc < 3){
    showUsage();
    exit(1);
  }
  options->kernel = 0;
  options->kernel_type = 0;
  options->sigma = 2.0;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->normalize = false;
  options->quietMode = false;
  options-> output_id = "";
  options->oneGramBySize = false;
  int i=3;
  while (i<argc){
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->kernel = atoi(argv[i+1]);
	  i+=2;
	}
      if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->kernel_type = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-N", 2) == 0)
	{
	  options->normalize = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-o", 2) == 0)
	{
	  options->output_id = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-gs", 3) == 0)
	{
	  options->oneGramBySize = true;
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
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  Options options;
  readOptions(argc, argv, &options);
  
  MoleculeGraph::initTable();
  //trainset : Base d'apprentissage
  MoleculesDataset trainset (options.trainset_path, options.dataset_file);
  trainset.computeLabeledSpectrums();
  trainset.computeRelevantCycles();
  trainset.computeCyclesSpectrums();


  double (*kernel)(double,double,double);
  switch(options.kernel_type){
  case 0:
    kernel = &intersectionKernel;
    break;
  case 1:
    kernel = &gaussianKernel;
    break;
  case 2:
    kernel = &innerKernel;
    break;
  case 3:
    kernel = &polynomialKernel;
    break;
  case 4:
    kernel = &gaussianKernel;
    break;
  default:
    cerr << "Undefined kernel."<< endl;
  }
  //vector< cimg_library::CImg<double> > Kks;
  vector< pair<string,int> > k_to_treelet;
  treelet_spectrum * distribution;
  
  int nb_treelets = 0;
  int N = trainset.size();
  char * spectrum_type;
  if(options.kernel == 0){
    distribution = trainset.getTreeletDistribution(always_true);
    spectrum_type = "labeled_spectrum";
  }  else if(options.kernel == 1){
    distribution = trainset.getCycleTreeletDistribution(always_true);
    spectrum_type = "cycle_spectrum";
  }   
  cout << "Nb Treelets : " << trainset.getNbTreelets() << endl;
  
  stringstream gram_filename;
  gram_filename << "gram_" << options.output_id << ".mat";
  ofstream outfile_gram (gram_filename.str().c_str());

  /*Parcours de tous les motifs de treelets*/
  if(!options.oneGramBySize)
    {
      for(int treelet_type=0;treelet_type < SIZE_SPECTRUM; treelet_type++){ 
	treelet_spectrum::iterator it = distribution[treelet_type].begin();
	for(;it != distribution[treelet_type].end();it++)
	  {
	    string treelet = it->first;
	    //cout <<  MoleculeGraph::translateTreeletCode(treelet) << " : " << it->second << endl;
	    CImg<double> Kk(N,N,1,1,0.0);
	    //XXX:Symétrique, à optimiser
	    /* Calcul sous matrice*/
	    for(int i=0;i<N;i++)
	      for(int j=0;j<N;j++)
		{
		  //XXX: We assume that the spectrum must be computed
		  pandore::Collection* c1 = trainset.getCollection(i);
		  pandore::Collection* c2 = trainset.getCollection(j);
		  treelet_spectrum ** spec_1 = 
		    (treelet_spectrum **) c1->GETVALUE(spectrum_type,Long);
		  treelet_spectrum ** spec_2 = 
		    (treelet_spectrum **) c2->GETVALUE(spectrum_type,Long);
		  
		  treelet_spectrum::iterator it1 = spec_1[treelet_type]->find(treelet);
		  treelet_spectrum::iterator it2 = spec_2[treelet_type]->find(treelet);
		  
		  double fi_k, fj_k;
		  //fi_k = 0 si pas de treelets dans la molécule c1
		  (it1 == spec_1[treelet_type]->end())?fi_k=0:fi_k=it1->second; 
		  (it2 == spec_2[treelet_type]->end())?fj_k=0:fj_k=it2->second;
		  double val;
		  if(options.kernel_type==1 && (fi_k==0 || fj_k==0))
		    {
		      val=0;
		    }
		  else
		    {
		      val=kernel(fi_k,fj_k,options.sigma);
		    }

		  Kk(i,j) = val;
		}
	    // Normalisation Matrice de Gram
	    
	    if(options.normalize)
	      for (unsigned int i=0; i<N; ++i)
		for (unsigned int j=0; j<N; ++j)
		  if ((Kk(i,i) != 0) && (Kk(j,j) != 0))
		    if(i!=j)
		      Kk(i,j) = Kk(i,j)/(sqrt(Kk(i,i)*Kk(j,j)));
	    
	    if(options.normalize)
	      for (unsigned int i=0; i<N; ++i)
		if (Kk(i,i) != 0)
		  Kk(i,i) = Kk(i,i)/(sqrt(Kk(i,i)*Kk(i,i)));
	    

	    
	    /* Verif semi defini positivité */
	    if(! options.quietMode){
	      CImg<double> eigvals;
	      CImg<double> eigvects;
	      Kk.symmetric_eigen(eigvals,eigvects);
	      cout << "Min eigenvalue of " << nb_treelets << " = " << eigvals[N-1] << endl;
	      double gram_sum = Kk.sum();
	      cout << "Sum of " << nb_treelets << " = " << gram_sum << endl;
	      /* Verif symétrie ok */
	      for (unsigned int i=0; i<N; ++i)
		for (unsigned int j=0; j<N; ++j)
		  if(Kk(i,j) != Kk(j,i))
		    cout << "Caramba : " << Kk(i,j) - Kk(j,i) << " " << Kk(i,j) << " " << Kk(j,i) <<  endl;
	      
	    }
	    /* Ecriture Fichier */
	    for (unsigned int i=0; i<N; ++i){
	      for (unsigned int j=0; j<N; ++j){
		outfile_gram << Kk(i,j);
		outfile_gram.flush();
		if(j!=N-1)
		  outfile_gram.flush() << ", ";
	      }
	      outfile_gram << endl;
	    }
	    // Kks.push_back (Kk);
	    k_to_treelet.push_back(pair<string,int>(treelet,treelet_type));
	    nb_treelets ++;
	    cout << nb_treelets << "/ " << trainset.getNbTreelets() << endl;
	  }
      }
    }
  else // If oneGramBySize
    {
      for(int treelet_type=0;treelet_type < 5; treelet_type++){
	treelet_spectrum::iterator it = distribution[treelet_type].begin();
	CImg<double> Kk(N,N,1,1,0.0);
	for(;it != distribution[treelet_type].end();it++)
	  {
	    string treelet = it->first;
	    /* Calcul sous matrice*/
	    for(int i=0;i<N;i++)
	      for(int j=i;j<N;j++)
		{
		  //XXX: We assume that the spectrum must be computed
		  pandore::Collection* c1 = trainset.getCollection(i);
		  pandore::Collection* c2 = trainset.getCollection(j);
		  treelet_spectrum ** spec_1 = 
		    (treelet_spectrum **) c1->GETVALUE(spectrum_type,Long);
		  treelet_spectrum ** spec_2 = 
		    (treelet_spectrum **) c2->GETVALUE(spectrum_type,Long);
		  
		  treelet_spectrum::iterator it1 = spec_1[treelet_type]->find(treelet);
		  treelet_spectrum::iterator it2 = spec_2[treelet_type]->find(treelet);
		  
		  double fi_k, fj_k;
		  //fi_k = 0 si pas de treelets dans la molécule c1
		  (it1 == spec_1[treelet_type]->end())?fi_k=0:fi_k=it1->second; 
		  (it2 == spec_2[treelet_type]->end())?fj_k=0:fj_k=it2->second;
		  double val;
		  if(options.kernel_type==1 && (fi_k==0 || fj_k==0))
		    {
		      val=0;
		    }
		  else
		    {
		      val=kernel(fi_k,fj_k,options.sigma);
		    }
		  Kk(i,j) += val;
		  if(i!=j)
		    Kk(j,i) += val;
		}
	  }
	// Normalisation Matrice de Gram
	if(options.normalize)
	  for (unsigned int i=0; i<N; ++i)
	    for (unsigned int j=0; j<N; ++j)
	      if ((Kk(i,i) != 0) && (Kk(j,j) != 0))
		if(i!=j)
		  Kk(i,j) = Kk(i,j)/(sqrt(Kk(i,i)*Kk(j,j)));
	
	if(options.normalize)
	  for (unsigned int i=0; i<N; ++i)
	     if (Kk(i,i) != 0)
	       Kk(i,i) = Kk(i,i)/(sqrt(Kk(i,i)*Kk(i,i)));

	    /* Ecriture Fichier */
	    for (unsigned int i=0; i<N; ++i){
	      for (unsigned int j=0; j<N; ++j){
		outfile_gram << Kk(i,j);
		outfile_gram.flush();
		if(j!=N-1)
		  outfile_gram.flush() << ", ";
	      }
	      outfile_gram << endl;
	    }
	    k_to_treelet.push_back(pair<string,int>("size:",treelet_type));
	    nb_treelets ++;
      }
    }
  outfile_gram.close();
  
  /* Ecriture autres fichiers */
 stringstream treelets_filename,
	treelets_cool_filename, y_filename;
  if(!options.oneGramBySize)
    {
      treelets_filename << "treelets_" << options.output_id << ".mat";
      ofstream outfile_treelets (treelets_filename.str().c_str());
      treelets_cool_filename << "treelets_cool_" << options.output_id << ".mat";
      ofstream outfile_treelets_cool (treelets_cool_filename.str().c_str());
      for(int k=0;k<nb_treelets;k++)
	{
	  outfile_treelets << k_to_treelet[k].second << "-" << k_to_treelet[k].first << endl;
	  if(options.kernel == 0){
	    outfile_treelets_cool << k_to_treelet[k].second << "-" << 
	      MoleculeGraph::translateTreeletCode(string(k_to_treelet[k].first)) << endl;
	  }else {
	    outfile_treelets_cool << k_to_treelet[k].second << "-" << 
	      MoleculeGraph::translateTreeletCycleCode(string(k_to_treelet[k].first)) << endl;
	  }
	}

      outfile_treelets.close();
      outfile_treelets_cool.close();
    }

  cout << "gram,treelets and treelets_cool files ok."<< endl;
  y_filename << "y_"<< options.output_id << ".mat" ;
  ofstream outfile_y (y_filename.str().c_str());
  for (int i=0; i<N; ++i)
      outfile_y << trainset.getParameter(i) <<  endl;
  outfile_y.close();
  cout << "y file ok."<< endl;
  
  return 0;
}


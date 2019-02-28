/*
 * @file moleculesClassification.cpp
 *
 * This program performs a SVM classification on molecules datasets.
 *
 *
 * @version 1.1.0 (2010-07-10)
 */

#include <iostream>

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
#include "GraphFileKernel.h"
using namespace std;
using namespace pandore;

//XXX: Code a factoriser
bool always_true(double i)
{
  return true;
}

map<char *, double,bool (*)(char *, char *)>* vit_to_weight(vector<char *> vit[SIZE_SPECTRUM])
{
  bool(*fn_pt)(char*,char*) = string_utils::keyComp;
  map<char *, double,bool (*)(char *, char *)> * weights = 
    new map<char *, double,bool (*)(char *, char *)>[SIZE_SPECTRUM];
  
  
  for(int k=0;k<SIZE_SPECTRUM;k++)
    {
      weights[k] = map<char *, double, bool (*)(char *, char *)> (fn_pt); 
      for(int i=0;i<vit[k].size();i++)
	weights[k].insert(pair<char*,double>(vit[k][i],1.0));
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
  char * file_gram;
  double c;
  bool quietMode;
  int useTestSet;
  bool normalize;
};

void showUsage ()
{
  cerr << "Usage : moleculesClassification trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-c C : Set the SVM parameter C (default 1.0)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-k : Gram Matrix File" << endl;
  cerr << "-t nb : Use nb last as test set, others as training set" << endl;
  cerr << "-r nb : Use 1/nb last as test set, others as training set" << endl;
  cerr << "-n : Normalize gram matrix" << endl;
  
}

void readOptions (int argc, char** argv, Options* options)
{
  if (argc < 3)
    {
      showUsage();
      exit(1);
    }
	
  // The default values are fixed
	
  options->c = 1.0;
  options->quietMode = false;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->file_gram = "";
  options->useTestSet = 0;
  options->normalize = 0;

  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->file_gram = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->useTestSet = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-n", 2) == 0)
	{
	  options->normalize = true;
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
  MoleculesDataset trainset (options.trainset_path, options.dataset_file);
  // We compute the Gram matrix of the training set
  GraphKernel* kgraph =  new GraphFileKernel(options.file_gram);
  trainset.computeGramMatrix(kgraph, options.normalize);
  trainset.showGramMatrixRaw();

  double fp=0,tp=0,fn=0,tn =0 ;
  double accuracy = 0.0;
  unsigned int N;

  double fpT=0,tpT=0,fnT=0,tnT =0 ;
  double accuracyT = 0.0;
  unsigned int NT;

  if(options.useTestSet ==0){
    // Leave-one-out cross-validation
    N = trainset.size();
    
    double  c;
    Collection* col;
    for (unsigned int i=0; i<N; ++i)
      {
	//On supprime le premier
	c = trainset.getParameter(0);
	col = trainset.getCollection(0);
	trainset.delete_first();
	
	SVM svm (trainset, kgraph, SVM::C_SVC, options.c);
	
	double x = svm.predict(col);
	
	if(!options.quietMode)
	  cout << "Molecule " << i << " : " << c << " --> " << x;
      
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
	
	else
	  if(!options.quietMode)
	    cout << "****" ;
      if(!options.quietMode)
	cout << endl;
      
      //On le rajoute a la fin
      trainset.add (col, c, kgraph);
      
      }
    
    if(!options.quietMode)
    cout << endl;
  }
  else 
    {
      N = options.useTestSet;
      double  *c = new double[N];
      Collection** col = new Collection*[N];
      //Nouveau dataset de train sur les 250 premieres
      for(unsigned int i=0;i<N;i++){
	c[i] = trainset.getParameter(trainset.size() -1);
	col[i] = trainset.getCollection(trainset.size() -1);
	trainset.delete_last();
      }
      SVM svm (trainset, kgraph, SVM::C_SVC, options.c);
      for(unsigned int i=0;i<N;i++){
	double x = svm.predict(col[i]);
	
	if(!options.quietMode)
	  cout << "Molecule " << i << " : " << c[i] << " --> " << x;
	else
	  if(!options.quietMode)
	    cout << "****" ;
	if(!options.quietMode)
	  cout << endl;
	
	if (c[i] == x)
	  ++accuracy;
	if ((c[i] == 1) && (x == 1))
	  tp ++;
	if ((c[i] == 1) && (x == 0))
	  fn ++;
	if ((c[i] == 0) && (x == 1))
	  fp ++;
	if ((c[i] == 0) && (x == 0))
	  tn ++;
      }

      NT=trainset.size();
      for(unsigned int i=0;i<NT;i++){
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
    }

  cout << "Accuracy Train : " << 100*accuracyT/NT << "% (" << accuracyT << " / " << NT << ")" << endl;
  double precisionT = tpT/(tpT+fpT);
  double recallT = tpT/(tpT+fnT);
  double f_scoreT = 2*(precisionT*recallT) / (precisionT+recallT);
  cout << precisionT << "," << recallT << "," << fpT << "," << tnT << ","<< tpT << "," << fnT << endl;
  cout << "F Score : " << f_scoreT << endl;


  cout << "Accuracy : " << 100*accuracy/N << "% (" << accuracy << " / " << N << ")" << endl;
  double precision = tp/(tp+fp);
  double recall = tp/(tp+fn);
  double f_score = 2*(precision*recall) / (precision+recall);
  cout << precision << "," << recall << "," << fp << "," << tn << ","<< tp << "," << fn << endl;
  cout << "F Score : " << f_score << endl;



  trainset.destroy();
  delete kgraph;
  
  return 0;
}

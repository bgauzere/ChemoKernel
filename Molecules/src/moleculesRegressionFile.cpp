/*
 * @file moleculesRegressionDirectory.cpp
 *
 * @description This program performs a Kernel Ridge regression on gram matrix files present in a directory.
 *
 * @author Benoit Gauzere <benoit.gauzere@ensicaen.fr>
 */

#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <pandore.h>

#include "CImg.h"
#include "MoleculesDataset.h"
#include "MoleculeGraph.h"
#include "GraphFileKernel.h"
#include "KernelRidge.h"
#include "SVM.h"

using namespace std;
using namespace pandore;

double predict (MoleculesDataset & trainset, GraphKernel* kgraph, double lambda, Collection* col,int methode,double c,double epsilon)
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
      delete kr;
    }
 
  return res;
}

void leaveOneOut (MoleculesDataset & dataset, GraphFileKernel* kgraph, double alpha, bool quietMode,int methode,double c,double epsilon)
{	
  // KRegression* kr = new KernelRidge (kgraph, &dataset, alpha);
  int N = dataset.size();

  //Stats data
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;
  double * bp_exps = new double[N];
  for (int i=0; i<N; ++i)
    {
      double bp_exp = dataset.getParameter(0);
      Collection* col = dataset.getCollection(0);
      char * filename = dataset.getMoleculeFilename(i);

      dataset.delete_first();
      double bp = predict (dataset, kgraph, alpha, col, methode, c, epsilon);
      double err = bp_exp - bp;
      err_quad += err * err;
      if(!quietMode){
	cout << "Molecule " << i << "(" << filename << ") : " << bp << " (" << bp_exp << ")";
	cout << "\t --> \t" << abs(err)  << "°" << endl;
      }
      err_moy += abs(bp_exp - bp);
      bp_mean += bp_exp; 
      bp_exps[i] = bp_exp;		
      dataset.add (col, bp_exp, kgraph);
    }
  
  bp_mean /= N;
  std_dev = sqrt(err_quad / N); 
  err_moy /= N;

  double r_den = 0.0;
  for(int i =0;i < N; i++){
    r_den += pow(bp_exps[i] - bp_mean,2);
  }
  double r = sqrt(1- (err_quad) / (r_den));
  cout << "\t Average error : " << err_moy  << endl;
  cout << "\t Standard Deviation : " << std_dev << endl;
  cout << "\t R : " << r << endl;
  delete [] bp_exps;
}

struct Options
{
  char* trainset_path;
  char* dataset_file;
  char * file_gram;
  double alpha;
  bool quietMode;
  int useTestSet;
  int useTestSetAlt;
  int regMethode;
  double c;
  bool normalize;
  double delta;
  double epsilon;
  int epsiPropSD;
};

void usage(char * s)
{
  cerr << "Usage : " << s << " trainset_path dataset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-k : path to gram matrix file" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-t nb : Use 1 molecule / nb as test, others as training set" << endl;
  cerr << "-s nb : Use nb first as train set, others as test set" << endl;
  cerr << "-reg regression : Choisis différentes méthodes de régression" << endl;
  cerr << "\t1 -- Kernel Ridge" << endl;
  cerr << "\t2 -- SVM" << endl;
  cerr << "-a alpha : Set the parameter of the regression for the Kernel Ridge (default 1.0)" << endl;
  cerr << "-c C SVM : Paramètre C du SVM" << endl; 
  cerr << "-n : Normalize the gram matrix" << endl; 
  cerr << "-d delta : Apply a RBF kernel on the gram matrix of parameter delta" <<endl;
  cerr << "-eps epsilon : Paramètre Epsilon du SVM " << endl;
  cerr << "-epsSD : Epsilon est prit comme proportionnel (à l'inverse de la valeur passer par -eps ) à l'écart type du jeu de données" << endl;

}

void readOptions (int argc, char** argv, Options* options)
{
  cout << endl;
  if (argc < 3)
    {
      usage(argv[0]);
      exit(1);
    }
	
  options->file_gram = "";
  options->alpha = 1.0;
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->quietMode = false;
  options->useTestSet = 0;
  options->useTestSetAlt = 0;
  options->regMethode=1;
  options->c=1.0;
  options->normalize=false;
  options->delta=0.0;
  options->epsilon=1;
  options->epsiPropSD=0;

  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->file_gram = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-a", 2) == 0)
	{
	  options->alpha = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->useTestSet = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->useTestSetAlt = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-reg", 4) == 0)
	{
	  options->regMethode = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-c", 2) == 0)
	{
	  options->c = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-n", 2) == 0)
	{
	  options->normalize = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->delta = atof(argv[i+1]);
	  if(options->delta!=0.0)
	    options->normalize = true;
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
  //Print command for archiving results
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }

  Options options;
  readOptions(argc, argv, &options);
  MoleculeGraph::initTable();
  MoleculesDataset * dataset = new MoleculesDataset(options.trainset_path, options.dataset_file);

  if(options.epsiPropSD)
    options.epsilon=dataset->standardDeviation()/options.epsilon;
  
  if(options.regMethode==2)
    cout<<"Epsilon : "<<options.epsilon<<endl;

  if(!options.quietMode)
    {
      cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
      cout << "Number of molecules : " << dataset->size() << endl;
      cout << "Alpha : " << options.alpha << endl;
      cout << endl;
    }
  
  cout << "file :" << options.file_gram << " : " << endl;
  GraphFileKernel * kgraph = new GraphFileKernel(options.file_gram);
  if(options.normalize)
    kgraph->normalize();
  if(options.delta!=0.0)
    kgraph->RBF(options.delta);

  dataset->computeGramMatrix(kgraph,false);

  cout << "Dataset Computed" << endl;
  cout << "Positive Definite  ?" << dataset->isGramMatrixPD() << endl;
  
  
  //Stats data
  int N = dataset->size();
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double * err_quad_by_struct = new double[SIZE_SPECTRUM];
  double bp_mean = 0.0;
  double * bp_exps = new double[N];
  int nb_molecules = 0;
  for(int i = 0; i < N; i ++){
    bp_exps[i] = dataset->getParameter(i);
  }
  
  if ((options.useTestSet == 0) && (options.useTestSetAlt == 0))
    leaveOneOut(*dataset, kgraph, options.alpha,options.quietMode,options.regMethode,options.c,options.epsilon);
  else if(options.useTestSetAlt > 0) {
    //Stats data
    int nb_app = options.useTestSetAlt;
    //Dataset split
    list<int> to_test;
    list<Collection *> to_test_Collections;
    list<double> to_test_Parameters;
    
    for(int n = nb_app; n < N; n++){
      to_test.push_back(n);
      to_test_Collections.push_back(dataset->getCollection(n));
      to_test_Parameters.push_back(dataset->getParameter(n));
    }
    dataset->eraseSome(to_test);
    //The Dataset is splitted
    double * coeffs;
    dataset->computeGramMatrix(kgraph, false);
    cout << "Gram Matrix Positive definite ? ! " << dataset->isGramMatrixPD() << endl;;
    //Test molecules prediction
    list<int>::iterator it = to_test.begin();
    list<Collection*>::iterator it_col = to_test_Collections.begin();    
    list<double>::iterator it_param = to_test_Parameters.begin();
    
    int index = 0;
    for(;it != to_test.end();++it, ++it_col, ++it_param){
      double bp = predict(*dataset, kgraph, options.alpha,*it_col,options.regMethode,options.c,options.epsilon);
      double bp_exp = *it_param;
      double err =  bp_exp - bp;
      err_quad +=  err * err;
      cout << bp << " (" << bp_exp << ")" "\t --> \t" << abs(err)  << "°" << endl;
      err_moy += abs(bp_exp - bp);
      bp_mean += bp_exp; 
      nb_molecules ++;
    } 
    bp_mean /= nb_molecules;
    std_dev = sqrt(err_quad / nb_molecules); 
    err_moy /= nb_molecules;
  
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
  }
  else {
    int step = options.useTestSet;
    for(int k = 0; k < step; k++) {
      cout << k << "/" << step<< endl;
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
      dataset->computeGramMatrix(kgraph, false);
      cout << "Gram Matrix Positive definite ? ! " << dataset->isGramMatrixPD() << endl;;
      //Test molecules prediction
      list<int>::iterator it = to_test.begin();
      list<Collection*>::iterator it_col = to_test_Collections.begin();    
      list<double>::iterator it_param = to_test_Parameters.begin();
      
      for(;it != to_test.end();++it, ++it_col, ++it_param){
	double bp = predict(*dataset, kgraph, options.alpha,*it_col,options.regMethode,options.c,options.epsilon);
	double bp_exp = *it_param;
	double err =  bp_exp - bp;
	err_quad +=  err * err;
	cout << bp << " (" << bp_exp << ")" "\t --> \t" << abs(err)  << "°" << endl;
	err_moy += abs(bp_exp - bp);
	bp_mean += bp_exp; 
	nb_molecules ++;
      }
      //Dataset cleaning
      delete dataset;
      dataset = new MoleculesDataset (options.trainset_path, options.dataset_file);
    }
    bp_mean /= nb_molecules;
    std_dev = sqrt(err_quad / nb_molecules); 
    err_moy /= nb_molecules;
  
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
  }
  
  delete kgraph;
  dataset->destroy();
  return 0;
}
      

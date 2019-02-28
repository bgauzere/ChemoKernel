/*
 * @file moleculesRegressionFileTrainAndTestSet.cpp
 *
 * @description This program performs a SVM regression on gram matrix files, split in a train test and a test set.
 *
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr
 */

#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <pandore.h>
#include <fstream>

#include "CImg.h"
#include "MoleculesDataset.h"
#include "MoleculeGraph.h"
#include "GraphFileKernel.h"
#include "KernelRidge.h"
#include "SVM.h"

using namespace std;
using namespace pandore;

double predict (MoleculesDataset & trainset, GraphKernel* kgraph, Collection* col,double c,double epsilon)
{ 
  double res;
  SVM svm (trainset, kgraph, SVM::EPSILON_SVR,c,0.5,epsilon);
  res=svm.predict(col);
 
  return res;
}

double leaveOneOut (MoleculesDataset & dataset, GraphFileKernel* kgraph, bool quietMode,double c,double epsilon)
{	
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
      // cout<<"K(0,0)="<<(*kgraph)(dataset.getCollection(0),dataset.getCollection(0))<<endl;
      // cout<<"K(1,0)="<<(*kgraph)(dataset.getCollection(1),dataset.getCollection(0))<<endl;
      // cout<<"K(2,0)="<<(*kgraph)(dataset.getCollection(2),dataset.getCollection(0))<<endl;
      // cout<<"K(3,0)="<<(*kgraph)(dataset.getCollection(3),dataset.getCollection(0))<<endl;
      // cout<<"K(i,0)="<<(*kgraph)(col,dataset.getCollection(0))<<endl;
      // cout<<"K(i,1)="<<(*kgraph)(col,dataset.getCollection(1))<<endl;


      // cout<<"----------------"<<i<<"--------------Gram:"<<endl;
      // dataset.showGramMatrixRaw();

      double bp = predict (dataset, kgraph, col, c, epsilon);
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
  double q = 1- (err_quad) / (r_den);

  cout << "Average Learning error : " << err_moy  << endl;
  cout << "Standard Deviation Learning : " << std_dev << endl;
  cout << "q² : " << q << endl;

  delete [] bp_exps;
  return q;
}

struct Options
{
  char* trainset_path;
  char* dataset_file;
  char* testset_file;
  char * file_gram;
  bool quietMode;
  double c;
  bool normalize;
  double delta;
  double epsilon;

};

void usage(char * s)
{
  cerr << "Usage : " << s << " trainset_path dataset_file testset_file [options]" << endl;
  cerr << "options:" << endl;
  cerr << "-k : path to gram matrix file (contain in the order the train and the test data)" << endl;
  cerr << "-q : Quiet mode. Displays only the percentage of good classification" << endl;
  cerr << "-c C SVM : Paramètre C du SVM" << endl; 
  cerr << "-n : Normalize the gram matrix" << endl; 
  cerr << "-d delta : Apply a RBF kernel on the gram matrix of parameter delta" <<endl;
  cerr << "-e epsilon : Paramètre Epsilon du SVM" <<endl;

}

void readOptions (int argc, char** argv, Options* options)
{
  cout << endl;
  if (argc < 4)
    {
      usage(argv[0]);
      exit(1);
    }
	
  options->file_gram = "";
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->testset_file = argv[3];
  options->quietMode = false;
  options->c=1.0;
  options->normalize=false;
  options->delta=0.0;
  options->epsilon=0.1;

  int i=4;
  while (i<argc)
    {
      if (strncmp(argv[i], "-k", 2) == 0)
	{
	  options->file_gram = argv[i+1];
	  i+=2;
	}
      else if (strncmp(argv[i], "-q", 2) == 0)
	{
	  options->quietMode = true;
	  i+=1;
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
      else if (strncmp(argv[i], "-e", 2) == 0)
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

  if(!options.quietMode)
    {
      cout << "Dataset : " << options.trainset_path << "/" << options.dataset_file << endl;
      cout << "Number of molecules : " << dataset->size() << endl;
      cout << endl;
    }
  
  cout << "file :" << options.file_gram << " : " << endl;
  GraphFileKernel * kgraph = new GraphFileKernel(options.file_gram);
  if(options.normalize)
    kgraph->normalize();
  if(options.delta!=0.0)
    kgraph->RBF(options.delta);


  // Recupération des noms des molécules du testSet
  vector<char *> testFile;
  char* testset;	
  if (options.trainset_path[strlen(options.trainset_path)-1] != '/')
    {
      testset = new char[strlen(options.trainset_path)+strlen(options.testset_file)+2];
      sprintf(testset, "%s/%s", options.trainset_path, options.testset_file);
    }
  else
    {
      testset = new char[strlen(options.trainset_path)+strlen(options.testset_file)+1];
      sprintf(testset, "%s%s", options.trainset_path, options.testset_file);
    }
  int sizeTest=0;
  ifstream f (testset);
  delete [] testset;
  long id = 0;
  if (f)
    {
      string s;
      while (getline(f, s))
	if (s[0] != '#')
	  {
	    const char * s_c_str = s.c_str();
	    const char * pch=strrchr(s_c_str,' ');
	    char* ctfile = new char [pch-s_c_str+1];
	    ctfile = strncpy(ctfile,s_c_str,pch-s_c_str);
	    ctfile[pch-s_c_str] = '\0';
	    testFile.push_back(ctfile);
	    sizeTest++;
	  }
    }
  f.close();

  // Recupération des données des molécules du testSet et suppresion du dataset
  double * testParameter= new double[sizeTest];
  Collection ** testCol = new Collection *[sizeTest];
  
  int current=0;
  for(vector<char *>::iterator it = testFile.begin();it != testFile.end(); ++it)
    {
      unsigned int num=dataset->getNumberByFilename(*it);
      testParameter[current] = dataset->getParameter(num);
      testCol[current] = dataset->getCollection(num);
      current++;
    }  
  dataset->eraseByFilename(testFile);


  dataset->computeGramMatrix(kgraph,false);

  // cout << "Dataset Computed" << endl;
  // cout << "Positive Definite  ?" << dataset->isGramMatrixPD() << endl;

  // cout<<"---------------Gram:"<<endl;
  // dataset->showGramMatrixRaw();



  // cout<<"Gram after destruction  of the test set:"<<endl;
  // dataset->showGramMatrixRaw();
  // cout << "Positive Definite  ?" << dataset->isGramMatrixPD() << endl;

  // Validation
  double q = leaveOneOut(*dataset, kgraph,options.quietMode,options.c,options.epsilon);
 
  // Test
  double err_moy = 0;
  double std_dev = 0.0;
  double err_quad = 0.0;
  double bp_mean = 0.0;
  double bp_exp_mean = 0.0;
  double * bp_exps = new double[sizeTest];
  double * bps = new double[sizeTest];
  for (int i=0; i<sizeTest; ++i)
    {
      double bp_exp = testParameter[i];
      Collection* col = testCol[i];
      char * filename = testFile[i];

      double bp = predict (*dataset, kgraph, col, options.c,options.epsilon);
      double err = bp_exp - bp;
      err_quad += err * err;
      if(!options.quietMode){
      	cout << "Test Molecule " << i << "(" << filename << ") : " << bp << " (" << bp_exp << ")";
      	cout << "\t --> \t" << abs(err)  << "°" << endl;
      }
      err_moy += abs(bp_exp - bp);
      bp_mean += bp_exp; 
      bp_exp_mean += bp_exp; 
      bp_exps[i] = bp_exp;
      bps[i]=bp;
    }
  
  bp_mean /= sizeTest;
  bp_exp_mean /= sizeTest;
  std_dev = sqrt(err_quad / sizeTest); 
  err_moy /= sizeTest;

  double r_den = 0.0;
  for(int i =0;i < sizeTest; i++){
    r_den += pow(bp_exps[i] - bp_mean,2);
  }

  double r_den2 = 0.0;
  for(int i =0;i < sizeTest; i++){
    r_den2 += pow(bps[i] - bp_mean,2);
  }
  double r_num =0.0;
  for(int i =0;i < sizeTest; i++){
    r_num += (bps[i] - bp_mean)*(bp_exps[i] - bp_exp_mean);
  }
  double r2 = r_num / (sqrt(r_den*r_den2));

  cout << "\t Average error : " << err_moy  << endl;
  cout << "\t Standard Deviation : " << std_dev << endl;

  cout << "R : " << r2 << endl;
  cout << "Norme : " << sqrt(q*q+r2*r2) << endl;

  delete [] bp_exps;
  delete [] bps;
  delete [] testParameter;
  delete [] testCol;

  testFile.clear();
  
  delete kgraph;
  dataset->destroy();
  return 0;
}

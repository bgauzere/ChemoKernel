/*
 * @file moleculesVectorFileClassification.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Jun 21 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
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
#include "VectorFileKernel.h"
#include "IntersectionKernel.h"
#include "GaussianKernel.h"
#include "InnerProductKernel.h"
using namespace std;
using namespace pandore;

//XXX: Code a factoriser
bool always_true(double i)
{
  return true;
}

struct Options
{
  char* dataset_path;
  char* dataset_file;
  char * vector_file;
  double c;
  bool quietMode;
  int useTestSet;
  int dim;
  int K;
  double sigma;
};


void usage (char * s)
{
  cerr << "Usage : "<< s <<  "dataset_path dataset_file" << endl;
  cerr << "options:" << endl;
  cerr << "-v vector_file" << endl;
  cerr << "-K kernel_type : Set the type of kernel between the two spectrums (Treelet Kernel)" << endl;
  cerr << "\t0 -- Intersection Kernel" << endl;
  cerr << "\t1 -- Gaussian Kernel" << endl;
  cerr << "\t2 -- Inner Product Kernel" << endl;
  cerr << "\t3 -- Polynomial Kernel" << endl;

  cerr << "-s sigma : Set the parameter for gaussian kernel or power for polynomial kernel" << endl;
  cerr << "-c C : SVM Regularization parameter" << endl;
  cerr << "-t nb : Use nb last as test set, others as training set" << endl;
  cerr << "-d dim : Vectors dimension" << endl;
  
}


void readOptions (int argc, char** argv, Options* options)
{
  if (argc < 3)
    {
      usage(argv[1]);
      exit(1);
    }
	
  // The default values are fixed
	
  options->c = 1.0;
  options->quietMode = false;
  options->dataset_path = argv[1];
  options->dataset_file = argv[2];
  options->vector_file = "";
  options->useTestSet = 0;
  options->sigma = 2.0;
  options->dim =1;

  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-v", 2) == 0)
	{
	  options->vector_file = argv[i+1];
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
      else if (strncmp(argv[i], "-K", 2) == 0)
	{
	  options->K = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->dim = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else
	{
	  cerr << "Unknown parameter : " << argv[i] << endl;
	  usage(argv[1]);
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
  MoleculesDataset dataset (options.dataset_path, options.dataset_file);
  cout << "Dataset Initialized "<< endl;
  Kernel *k=NULL;
  switch(options.K){
  case Kernel::IntersectionKernelType:
    k=new IntersectionKernel();
    break;
  case Kernel::GaussianKernelType:
    k=new GaussianKernel(options.sigma);
    break;
  case Kernel::InnerProductKernelType:
    k=new InnerProductKernel();
    break;
  default:
    cerr << "Kernel undefined" << endl;
    exit(EXIT_FAILURE);
  }

  GraphKernel* kgraph =  new VectorFileKernel(options.vector_file,&dataset,options.dim,k);
  cout << "Graph Kernel ok" << endl;
  dataset.computeGramMatrix(kgraph, false);
  // trainset.showGramMatrixRaw();

  double fp,tp,fn,tn =0 ;
  double accuracy = 0.0;
  unsigned int N;
  if(options.useTestSet ==0){
    // Leave-one-out cross-validation
    N = dataset.size();
    
    double  c;
    Collection* col;
    for (unsigned int i=0; i<N; ++i)
      {
	//On supprime le premier
	c = dataset.getParameter(0);
	col = dataset.getCollection(0);
	dataset.delete_first();
	
	SVM svm (dataset, kgraph, SVM::C_SVC, options.c);
	
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
	
      if(!options.quietMode)
	cout << endl;
      
      //On le rajoute a la fin
      dataset.add (col, c, kgraph);
      
      }
    
    if(!options.quietMode)
    cout << endl;
  }
  else 
    {
      N = options.useTestSet;
      double  *c = new double[N];
      Collection** col = new Collection*[N];;
      //Nouveau dataset de train sur les 250 premieres
      for(unsigned int i=0;i<N;i++){
	c[i] = dataset.getParameter(dataset.size() -1);
	col[i] = dataset.getCollection(dataset.size() -1);
	dataset.delete_last();
      }
      SVM svm (dataset, kgraph, SVM::C_SVC, options.c);
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
    }
  cout << "Accuracy : " << 100*accuracy/N << "% (" << accuracy << " / " << N << ")" << endl;

  dataset.destroy();
  delete kgraph;
  
  return 0;
}

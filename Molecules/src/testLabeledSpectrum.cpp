/*
 * @file testLabeledSpectrum.cpp
 * @author Benoit Gauzere <<benoit.gauzere@ensicaen.fr>> 
 * @version     0.0.1 - Thu Feb  3 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * This programs the enumeration on Labeled treelet in a graph
 * All necessary references.
 *
 */

#include <iostream>
#include <pandore.h>
#include "MoleculesDataset.h"
#include "TreeletKernel.h"
#include "CImg.h"
#include "MoleculeGraphEditDistance.h"
#include "KEditDistance.h"
#include "InterTreeletKernel.h"
#include "LaplacianKernelOriginal.h"
#include "GaussianKernel.h"
using namespace std;
using namespace cimg_library;
using namespace pandore;

void usage (char * s)
{
  cerr << "Usage : "<< s << " path_to_dataset dataset_file" << endl;
  cerr << "options:" << endl;
}
bool always_true(double i)
{
  return true;
}

int main (int argc, char** argv)
{	
  MoleculeGraph::initTable();
  MoleculesDataset dataset (argv[1],argv[2]);
  dataset.computeLabeledSpectrums();
					      
  treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
  int nb_treelets =  dataset.getNbTreelets();
  cout << nb_treelets<< endl;
  dataset.computeSpecialVectors(&distrib, nb_treelets);
  cout << "Special Vectors computed " << endl;
  // for(int i=0;i<dataset.size();i++){
  //   CImg<double> * vec = dataset.getSpecialVector(i);
  //   for(int j=0;j<nb_treelets;j++)
  //     cout << (*vec)(j) << ",";
  //   cout << endl;
  // }
  // MoleculesDataset dataset (argv[1],argv[2]);
  // dataset.computeLabeledSpectrums();
  // dataset.printLabeledSpectrums();
  GraphEditDistance* edit = new MoleculeGraphEditDistance;
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
  cout << atof(argv[3]) << endl;
  //GraphKernel * kinter_treelets = new KEditDistance (edit, atof(argv[3]));
  GraphKernel * kinter_treelets = new LaplacianKernelOriginal (edit, dataset_treelets,atof(argv[3]),1,1);
  cout << "Laplacian Computed" << endl;
  Kernel * kinter_graphs = new GaussianKernel(1.0);
  cout << "TreeletKernel Initialized" << endl;
  GraphKernel * kgraph = new InterTreeletKernel(kinter_treelets,kinter_graphs,
						// &distrib,
						col_treelets,
						nb_treelets);
  CImg<double> * sim_treelets = ((InterTreeletKernel *)kgraph)->getSimTreelets();
  cout << "Gram inter treelets : "<< endl;
  for(int i=0;i<nb_treelets;i++){
    for(int j=0;j<nb_treelets;j++)
      cout << (*sim_treelets)(i,j)<< " ";
    cout << endl;
  }
  CImgDisplay main_disp(*sim_treelets,"Sim Treelets");
  cout << "Kernel Initialized " << endl; 
  dataset.computeGramMatrix(kgraph);
  cout << "Gram Matrix computed " << endl;
  cout << "Is Semi Positive Definite ?" << endl;
  cout << dataset.isGramMatrixPD() << endl;
  cout << "Sdp computed " << endl;
  dataset.showGramMatrixMatlab("gram_inter_treelets.mat");
  cout << "Done " << endl;
  return 0;
}

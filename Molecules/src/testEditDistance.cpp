/*
 * @file testEditDistance.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon May 21 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include "CImg.h"
#include "MoleculesDataset.h"
#include "MoleculeGraphEditDistanceMCS.h"
#include "TreeletEditDistance.h"
#include <pandore.h>

using namespace std;
using namespace cimg_library;
using namespace pandore;

void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

bool always_true(double i)
{
  return true;
}


int main (int argc, char** argv)
{	
  TreeletEditDistance::Init();
  MoleculeGraph::initTable();
  cout << TreeletEditDistance::maxPermsToCompute() << endl;
  exit(EXIT_SUCCESS);

  MoleculesDataset dataset (argv[1], argv[2]);
  dataset.computeLabeledSpectrums();
  treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
  int nb_treelets =  dataset.getNbTreelets();
  cout << nb_treelets<< endl;
  dataset.computeSpecialVectors(&distrib, nb_treelets);
  cout << "Special Vectors computed " << endl;

  double cout_substitution = 1;
  double cout_structurel = 6;

  MoleculeGraphEditDistanceMCS edit(cout_substitution,cout_structurel);
  Collection ** col_treelets = new Collection* [nb_treelets];
  string * codes_treelets = new string[nb_treelets];
  int * types_treelets = new int[nb_treelets];
  int cur_treelet = 0;
  for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
    {
      treelet_spectrum::iterator it = distrib[i].begin();
      for(;it != distrib[i].end();it++,cur_treelet++){
	//cout << i << "," << it->first << endl;
	Collection * tmp = TreeletEnumerator::TreeletToCollection(i,it->first);
	col_treelets[cur_treelet] = tmp;
	codes_treelets[cur_treelet] = it->first;
	types_treelets[cur_treelet] = i;
      }
    }
  TreeletEditDistance::Init();
  TreeletEditDistance t(cout_substitution,cout_structurel);
  CImg<double> ed(nb_treelets,nb_treelets);
  CImg<double> suboptimal_ed(nb_treelets,nb_treelets);
  for(int i=0;i<nb_treelets;i++)
    for(int j=i;j<nb_treelets;j++){
      ed(i,j) = ed(j,i) = t(types_treelets[i],codes_treelets[i],
		  types_treelets[j],codes_treelets[j]);
      double ed1 = edit(col_treelets[i],col_treelets[j]);
      double ed2 = edit(col_treelets[j],col_treelets[i]);
      suboptimal_ed(i,j) = suboptimal_ed(j,i) = (ed1<ed2)?ed1:ed2;
      if(suboptimal_ed(i,j)<ed(i,j)){
	cout << "(" <<i<<","<<j<<")"<<types_treelets[i]<< "-" << 
	  MoleculeGraph::translateTreeletCode(codes_treelets[i]);
	cout << " => " << types_treelets[j]<< "-" << MoleculeGraph::translateTreeletCode(codes_treelets[j]) << endl;
	cout << "ed=" << ed(i,j) << endl;
	cout << "sub optimal=" << suboptimal_ed(i,j) << endl;
      }
    }
  
  for(int i=0;i<nb_treelets;i++)
    for(int j=0;j<nb_treelets;j++)
      if(suboptimal_ed(i,j) !=  suboptimal_ed(j,i))
	{
	  cout << "(" <<i<<","<<j<<")"<<types_treelets[i]<< "-" << 
	    MoleculeGraph::translateTreeletCode(codes_treelets[i]);
	  cout << " => " << types_treelets[j]<< "-" << MoleculeGraph::translateTreeletCode(codes_treelets[j]) << endl;
	  cout << "sub optimal(i,j)=" << suboptimal_ed(i,j) << endl;
	  cout << "sub optimal(j,i)=" << suboptimal_ed(j,i) << endl;
	}
      
  ed.display();
  suboptimal_ed.display();
  CImg<double> diff(nb_treelets,nb_treelets);
  diff = suboptimal_ed - ed;
  //for(int i=0;i<nb_treelets;i++)
  //for(int j=0;j<nb_treelets;j++)
      //     diff(i,j) = (diff(i,j)<0)?-diff(i,j):0;
  diff.display();

  CImg<double> eigvals;
  CImg<double> eigvects;
  CImg<double> k(nb_treelets,nb_treelets);
  for(int i=0;i<nb_treelets;i++)
    for(int j=0;j<nb_treelets;j++)
      k(i,j) = exp(-pow(ed(i,j),2)/(5));
  k.display();
  k.symmetric_eigen(eigvals,eigvects);
  for (unsigned int i=0; i<nb_treelets; ++i)
    cout << eigvals[i] << endl;

  cout << "Sub optimal" << endl;
  for(int i=0;i<nb_treelets;i++)
    for(int j=0;j<nb_treelets;j++)
      k(i,j) = exp(-pow(suboptimal_ed(i,j),2)/(5));
  k.display();
  k.symmetric_eigen(eigvals,eigvects);
  for (unsigned int i=0; i<nb_treelets; ++i)
    cout << eigvals[i] << endl;

    

  return 0;
}

/*
 * @file moleculesCalculNbCycles.cpp
 *
 * @description This program computes the size of a Molecule for each Molecule in a Dataset.
 *
 * Usage : moleculesCalculNbCycles trainset_path dataset_file
 * 
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 * @version 1.0 (2010-10-25)
 */

#include <iostream>
#include <pandore.h>

#include "GraphKernel.h"
#include "KEditDistance.h"
#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "MoleculeGraph.h"
#include "utils.h"
#include "TreeletEnumeratorAugmentedCycles.h"
#include "separator.h"
using namespace std;
using namespace pandore;

void usage(char * s)
{
  cerr << "This program computes the number of rings of each molecule in the given dataset." << endl;
  cerr << "Usage : " << s << " dataset_path dataset_file" << endl;
}

#define NB_ARGS 2
int main (int argc, char** argv)
{	
#ifdef DEBUG
  cout << "toto" << endl;
#endif
  MoleculeGraph::initTable();
  if(argc != NB_ARGS + 1)
    {
      usage(argv[0]);
      return -1;
    }
  MoleculesDataset dataset (argv[1],argv[2]);
  int N = dataset.size();
  dataset.computeCCHypergraph();
  dataset.computeAugmentedCycles();
  dataset.computeAugmentedCycleSpectrums();
	
  for (int i=0; i<N; ++i)
    {
      Collection* col = dataset.getCollection(i);
      int nb_vertex = col->GETPARRAYSIZE("nodes", Collection);
      int nb_edges = col->GETARRAYSIZE("edges", Char);
      cout << dataset.getMoleculeFilename(i) << " ";
      cout <<  "Nb cyclomatique:" << nb_edges - nb_vertex + 1 << endl;
      
      treelet_spectrum ** spec;
      spec = dataset.getAugmentedCycleSpectrum(0);
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  {
	    treelet_spectrum::iterator it = spec[j]->begin();
	    for(;it != spec[j]->end();it ++)
	      {
		cout << j << ":" ;
		string code = it->first;
		string::iterator it_string;
		for(it_string = code.begin();it_string != code.end();it_string ++){
		  char c = *it_string;
		  if (c == ANGLE_CHAR_SEP)
		    cout << c;
		  else
		    cout << (unsigned int)(*it_string) << ",";
		}
		cout << " :" << it->second << endl;
	      }
	  }
    }
  
  return 0;
}

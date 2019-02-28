/*
 * @file moleculesCalculNbAtoms.cpp
 *
 * @description This program computes the size of a Molecule for each Molecule in a Dataset.
 *
 * Usage : moleculesCalculNbAtoms trainset_path dataset_file
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

using namespace std;
using namespace pandore;

void usage(char * s)
{
  cerr << "This program computes the size of each molecule in the given dataset." << endl;
  cerr << "Usage : " << s << " dataset_path dataset_file" << endl;
}

#define NB_ARGS 2
int main (int argc, char** argv)
{	
  if(argc != NB_ARGS + 1)
    {
      usage(argv[0]);
      return -1;
    }
  MoleculesDataset dataset (argv[1],argv[2]);
  int N = dataset.size();
	
  for (int i=0; i<N; ++i)
    {
      Collection* col = dataset.getCollection(i);
      cout << col->GETPARRAYSIZE("nodes", Collection) << endl;

    }
	
  return 0;
}

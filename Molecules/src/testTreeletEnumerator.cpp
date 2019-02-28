/*
 * @file testTreeletEnumerator.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Feb  7 2012
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
#include "MoleculesDataset.h"
#include "HorvathKernel.h"

using namespace std;
bool always_true(double i)
{
  return true;
}

void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{
  MoleculeGraph::initTable();
  MoleculesDataset dataset (argv[1], argv[2]);
  dataset.computeLabeledSpectrums();
  treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
  cout << "Treelets :" << endl;
  for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      treelet_spectrum::iterator it = distrib[treelet_type].begin();
      for(;it != distrib[treelet_type].end();it ++)
  	cout << treelet_type << "-" << MoleculeGraph::translateTreeletCode(it->first) << " : " << it->second<< endl;
    }
  delete [] distrib;
   
  return 0;
}

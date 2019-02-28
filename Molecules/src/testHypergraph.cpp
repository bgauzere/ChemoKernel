/*
 * @file testHypergraph.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Dec  5 2012
 * 
 * @todo  
 * @bug 
 *  
 * Test program for contracted cycles hypergraphs building
 * 
 *
 */

#include <iostream>
#include "MoleculeGraph.h"
#include "MoleculesDataset.h"

using namespace std;

void usage (char * s)
{
  cerr << "Usage : "<<s << " path_to_dataset dataset_file" << endl;
}

bool always_true(double i)
{
  return true;
}

int main (int argc, char** argv)
{	
  MoleculeGraph::initTable();
  MoleculesDataset dataset (argv[1],argv[2]);
  dataset.computeCCHypergraph();
  dataset.computeContractedCycleSpectrums();
  // dataset.computeSpectrums();
  // for (int i =0;i<dataset.size();i++){
  //   cout << "*******************************" << endl;
  //   cout << "\t Molecule "  << i << endl;
  //   cout << "*******************************" << endl;
  //   MoleculeGraph * mol = dataset.getMoleculeGraph(i);
  //   mol->computeContractedCycleGraphs();
  // }
  treelet_spectrum * distrib = dataset.getContractedCycleTreeletDistribution(&always_true);
  cout << "Treelets :" << endl;
  int nb_treelets = 0;
  int * histo = new int[10000];
  memset(histo, 0, 10000*sizeof(int));
  double max = 0;
  for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      treelet_spectrum::iterator it = distrib[treelet_type].begin();
      for(;it != distrib[treelet_type].end();it ++){
  	//cout << treelet_type << "-" << MoleculeGraph::translateTreeletCode(it->first) << " : " << it->second<< endl;
  	int cast_nb_occ = static_cast<int>(it->second);
  	histo[cast_nb_occ] ++;
  	max=(max>cast_nb_occ)?max:cast_nb_occ;
      }
    }
  cout << "Nb treelets = " << nb_treelets<< endl;
  for(int i=0;i<max+1;i++)
    cerr << i << " " << histo[i] << endl;
  delete [] distrib;
  cout << max << endl;
  return 0;
}

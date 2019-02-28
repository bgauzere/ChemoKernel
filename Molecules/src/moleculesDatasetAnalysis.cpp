/*
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
#include <vector>
#include "MoleculeGraph.h"
#include "MoleculesDataset.h"
#include "CycleKernel.h"

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
  int max_atoms = 0;
  int min_atoms = INT_MAX;
  int sum_atoms = 0;
  int sum_edges = 0;
  int sum_nbStereocenter=0;
  int max_sizeMSSG =0;
  int min_sizeMSSG =INT_MAX;
  int sum_sizeMSSG =0;

  for (unsigned int i =0;i<dataset.size();i++){
    MoleculeGraph * m = dataset.getMoleculeGraph(i);
    int m_atoms = m->nbAtoms();
    int m_edges = m->nbEdges();
    max_atoms = (m_atoms>max_atoms)?m_atoms:max_atoms;
    min_atoms = (m_atoms<min_atoms)?m_atoms:min_atoms;
    sum_edges = sum_edges + m_edges;
    sum_atoms = sum_atoms + m_atoms;

  }


  cout << "Nb mol : " << dataset.size() << endl;
  cout << "Ecart Type de la propriété : " << dataset.standardDeviation(true) << endl;
  cout << "Taille Moyenne : " << sum_atoms/(double)dataset.size() << endl;
  cout << "Degré Moyen : " << 2*sum_edges/(double)sum_atoms << endl;
  cout << "Max atomes : " << max_atoms << endl;
  cout << "Min atomes : " << min_atoms << endl;
  

  // dataset.computeLabeledSpectrums();


  // treelet_spectrum * distrib = dataset.getTreeletDistribution(&always_true);
  // cout << "Treelets :" << endl;
  // for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
  //   {
  //     treelet_spectrum::iterator it = distrib[treelet_type].begin();
  //     for(;it != distrib[treelet_type].end();it ++)
  // 	cout << treelet_type << "-" << MoleculeGraph::translateTreeletCode(it->first) << " : " << it->second<< endl;
  //   }
  // //Affichage molécule par molécule
  
  // treelet_spectrum ** distrib_by_treelet;// = new treelet_spectrum**[dataset.size()];
  // int N = dataset.size();
  // for(int i=0;i<N;i++){
  //   distrib_by_treelet = dataset.getTreeletSpectrum(i);
  //   cout << dataset.getMoleculeFilename(i) << endl;
  //   for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
  //     {
  // 	treelet_spectrum::iterator it = distrib_by_treelet[treelet_type]->begin();
  // 	for(;it != distrib_by_treelet[treelet_type]->end();it ++)
  // 	  cout << treelet_type << "-" << MoleculeGraph::translateTreeletCode(it->first) << " : " << it->second << endl;
  // 	//cout << "G" << treelet_type << " : " << it->second << endl;
  //     }
  // }
  
  
  // // cout << "Cycles :" << endl;
  // dataset.computeRelevantCycles();
  // dataset.computeCyclesSpectrums();
  // treelet_spectrum * distrib = dataset.getCycleTreeletDistribution(&always_true);
  // for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
  //   {
  //     treelet_spectrum::iterator it = distrib[treelet_type].begin();
  //     for(;it != distrib[treelet_type].end();it ++)
  // 	{//Parcours de tous les treelets:
  // 	  cout << treelet_type << "-" << 
  // 	    MoleculeGraph::translateTreeletCycleCode(it->first) << " : " << it->second << endl;
  // 	}
  //   }
  
  // treelet_spectrum *** distrib_by_treelet = new treelet_spectrum**[dataset.size()];
  // int N = dataset.size();
  // for(int i=0;i<N;i++){
  //   distrib_by_treelet[i] = dataset.getCycleTreeletSpectrum(i);
  // }
  // //   for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
  // //     {
  // // 	treelet_spectrum::iterator it = distrib_by_treelet[i][treelet_type]->begin();
  // // 	for(;it != distrib_by_treelet[i][treelet_type]->end();it ++)
  // // 	  cout << treelet_type << "-" << MoleculeGraph::translateTreeletCycleCode(it->first) << " : " << it->second << endl;
  // //     }
  // // }
  
  // for(int treelet_type =0;treelet_type<SIZE_SPECTRUM;treelet_type++)
  //   {
  //     treelet_spectrum::iterator it = distrib[treelet_type].begin();
  //     for(;it != distrib[treelet_type].end();it ++)
  // 	{//Parcours de tous les treelets:
  // 	  cout << treelet_type << "-" << it->first << " : " << it->second << "," ;
  // 	  cout << it->second/N << ", ";
  // 	  int nb_founds = 0;
  // 	  for(int i=0;i<N;i++){
  // 	    if((distrib_by_treelet[i][treelet_type])->find(it->first) != (distrib_by_treelet[i][treelet_type])->end())
  // 	      {
  // 		// cout << (distrib_by_treelet[i][treelet_type])->find(it->first)->second << ",";
  // 		nb_founds ++;
  // 	      }
  // 	  }
  // 	  cout << nb_founds  << endl;
  // 	}
  //   }
  // delete [] distrib;
  // // cout << "Simple Cycles" << endl;
  // // dataset.computeSimpleCycles(3);
  // // vector<string> codes = dataset.getSimpleCycleCodes(0);
  // // HorvathKernel * kgraph = new HorvathKernel();
  // // dataset.computeGramMatrix(kgraph,false);
  // // // dataset.showGramMatrixRaw();
  // // // for(int i=0;i<codes.size();i++)
  // // //   cout << MoleculeGraph::translateTreeletCycleCode(codes[i]) << endl;
  
  return 0;
}

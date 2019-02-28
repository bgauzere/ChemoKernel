/*
 * @file MoleculeSpectrumKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Apr  5 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "MoleculeSpectrumKernel.h"
#include <cassert>
#include <string>
using namespace std;

bool MoleculeSpectrumKernel::isAdmitted(string treelet, int treelet_type){
  for(unsigned int i = 0;i< _vit_list[treelet_type].size();i++)
    if(! _vit_list[treelet_type][i].compare(treelet))
      return true;
  return false;
}

void MoleculeSpectrumKernel::selectTreelets(std::vector<string> treelets[SIZE_SPECTRUM])
{
  for(int i=0;i<SIZE_SPECTRUM;i++)
      _vit_list[i] = vector<string>(treelets[i]);
  
}

double MoleculeSpectrumKernel::getWeight(string treelet, int treelet_type){
  map<string, double>::iterator elt = _weights_treelet[treelet_type].find(treelet);
  if(elt != _weights_treelet[treelet_type].end())
    return elt->second;
  else
    {
      cerr << "The weight for treelet G" << treelet_type << "-" << treelet 
	   << " is not defined." << endl;
      assert(0);
    }
}
 
void MoleculeSpectrumKernel::weightTreelets(treelet_spectrum treelets[SIZE_SPECTRUM])
{
  for(int i=0;i<SIZE_SPECTRUM;i++)
    _weights_treelet[i] = treelet_spectrum(treelets[i]);

}

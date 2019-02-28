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

#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>
#include "CImg.h"
#include "MoleculeGraph.h"
#include "MoleculesDataset.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculeGraphEditDistance.h"

using namespace std;
using namespace cimg_library;

void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{
  MoleculeGraph::initTable();
  MoleculesDataset dataset (argv[1], argv[2]);
  int n = dataset.size();
  MoleculeGraphEditDistance mcs;
  CImg<double> gram(n,n);
  for (int i = 0; i<n;i++)
    for (int j = i; j<n;j++)
      gram(i,j) =  gram(j,i) = mcs(dataset.getCollection(i), dataset.getCollection(j));
  
  ofstream outfile_gram (argv[3]);

  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int j=0; j<n; ++j)
	{
	  outfile_gram.precision(DBL_DIG);
	  outfile_gram <<  gram(i,j);
	  outfile_gram.flush();
	  if(j!=n-1)
	    outfile_gram.flush() << ", ";
	}
      outfile_gram << endl;
    }
  outfile_gram.close();

  return 0;
}

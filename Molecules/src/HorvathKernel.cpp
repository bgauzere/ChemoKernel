/*
 * @file HorvathKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Mar 19 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <vector>
#include <algorithm>

#include "HorvathKernel.h"
#include "MoleculeGraph.h"
using namespace pandore;
using namespace std;

double HorvathKernel::operator() (Collection* c1, Collection* c2)
{

  vector<string>  v_1 = (*(vector<string>*) c1->GETVALUE("simple_cycle_codes",Long));
  vector<string>  v_2 = (*(vector<string>*) c2->GETVALUE("simple_cycle_codes",Long));
  sort(v_1.begin(), v_1.end());
  sort(v_2.begin(), v_2.end());
  vector<string> result(v_1.size()+v_2.size());
  vector<string>::iterator it=set_intersection (v_1.begin(), v_1.end(),
						v_2.begin(), v_2.end(), result.begin());
  
  // cout << "intersection has " << int(it - result.begin()) << " elements.\n";
  return double(it - result.begin());

}


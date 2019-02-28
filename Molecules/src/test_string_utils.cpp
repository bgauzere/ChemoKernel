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

#include <cstdlib>
#include <iostream>
#include <cstring>
#include "string_utils.h"
#include "utils.h"
using namespace std;

void usage (char * s)
{
  cerr << "Usage : " << s << " chaine k " << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{	
  if (argc != 3)
    {
      usage(argv[0]);
      return -1;
    }

  char * seq = argv[1];
  int n = strlen(seq);
  int k = atoi(argv[2]);
  char ** perm_list;
  int ** perm_pos_init;
    
  char null_char = 0;
  int nb = utils::get_perm(seq, n, perm_list, perm_pos_init, k, null_char);

  cout << seq << endl;
  cout << "Nb perm : " << nb << endl;
  for(int i = 0; i < nb; i ++)
    {
      cout << i << " : ";
      for(int ii=0;ii<k;ii++)
	cout << perm_list[i][ii];
      cout << " [ ";
      for(int j=0;j<k;j++)
	cout << perm_pos_init[i][j] << " ";
      cout << "]" << endl;
    }
  for(int i=0;i<nb;i++){
    delete [] perm_list[i];
    delete [] perm_pos_init[i];
  }
  delete [] perm_pos_init;
  delete [] perm_list;
  return 0;
}

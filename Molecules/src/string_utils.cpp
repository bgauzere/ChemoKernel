/*
 * @file string_utils.cpp
 * @author Benoit Gauzere <<benoit.gauzere@ensicaen.fr>> 
 * @version     0.0.1 - Fri Feb  4 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include <cstdlib>
#include <cstring>
#include "utils.h"
#include "string_utils.h"

char * string_utils::reverse_string(const char * orig){
  int n = strlen(orig);
  char * reverse = new char[n+1];
  for(int i = 0; i<n; i++)
    reverse[i] = orig[n-i-1];
  reverse[n] = '\0';
    return reverse;  
}

bool string_utils::keyComp(char * l, char * r){
  return (strcmp(l,r) < 0);
}

bool string_utils::keyStringComp(std::string l, std::string r){
  return l.compare(r) < 0;
}

int string_utils::genericStringComp(const void * l, const void * r){
  return strcmp((const char *)l,(const char *)r);
}


char * string_utils::stringSort(const char * to_sort){
  char * sort_string = new char[strlen(to_sort)+1];
  sort_string = strcpy(sort_string,to_sort);
  qsort (sort_string,strlen(to_sort),sizeof(char),genericStringComp);
  return sort_string;  
}

int string_utils::get_perm(char * s, int seq_size, char *** perm_list, int *** perm_pos_init, int k)
{
  std::string seq = std::string(s);
  int nb_perm =  gsl_sf_choose(seq_size, k);
  /* Mem Alloc */
  (*perm_list) = new char*[nb_perm];
  (*perm_pos_init) = new int*[nb_perm];
  for(int i=0;i<nb_perm;i++)
    {
      (*perm_pos_init)[i] = new int[k];
      (*perm_list)[i] = new char[k+1];// \0
    }
  //Cas particulier n = k
  if (nb_perm ==1)
    {
      (*perm_list)[0] = strncpy((*perm_list)[0],s,k);
      (*perm_list)[0][k] = '\0';
      for(int i=0;i<k;i++)
	(*perm_pos_init)[0][i] = i;
    }
  else
    {
      int perm_pos = 0;
      string_utils::get_perm_rec(seq,seq_size,(*perm_list),(*perm_pos_init),k,&perm_pos,0,0);
    }
  return nb_perm;

}


void string_utils::get_perm_rec(std::string seq, int n, char ** perm_list, int ** perm_pos_init, 
				int k, int * perm_pos, int start_perm, int nb_mask)
{
  for(int i=start_perm; i<n; i++)
    {
      std::string seq_copy = std::string(seq);
      seq_copy[i] = '0';
      if ((n-nb_mask-1) == k)     //Permutation #perm_pos get !
	{
	  int pos=0;
	  for(int j=0; j < seq_copy.length(); j++) 
	    if(seq_copy[j] != '0')
	      {
		//j : position dans la chaine de caractere en entrée
		perm_pos_init[(*perm_pos)][pos] = j;
		perm_list[(*perm_pos)][pos] = seq_copy[j];
		pos ++;
	      }
	  perm_list[(*perm_pos)][pos] = '\0';
	  (*perm_pos) ++; //Permutation computed !
	}
      else if(n-nb_mask > k) //Pruning
	{
	  //seq_copy = strcpy(seq_copy, seq);
	  get_perm_rec(seq_copy,n, perm_list,perm_pos_init,k,perm_pos,i+1,nb_mask+1);
	}
      //delete [] seq_copy;
    }
}

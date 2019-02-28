
/**
 * @file utils.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Jun 13 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * Gnu Scientific Library required
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>
#include <iostream>
#include "CImg.h"
#include "pandore.h"
using namespace cimg_library;

static int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

static  double gsl_sf_choose (unsigned int n, unsigned int m){
  return factorial(n)/(factorial(m)*factorial(n-m));
}

class utils
{
private:
  template<typename T>
  static void get_perm_rec(const T * seq, int n, T ** perm_list, int ** perm_pos_init, 
			   int k, int * perm_pos, int start_perm, int nb_mask, int seq_size, T& null_element){
    for(int i=start_perm; i<n; i++)
      {
	T* seq_copy = new T[seq_size];
	for(int ii=0;ii<seq_size;ii++)
	  seq_copy[ii] = seq[ii];
	//memcpy(seq_copy,seq,seq_size*sizeof(T));
	seq_copy[i] = null_element;
	if ((n-nb_mask-1) == k)     //Permutation #perm_pos get !
	  {
	    int pos=0;
	    for(int j=0; j < seq_size; j++) 
	      if(seq_copy[j] != null_element)
		{
		  //j : position dans la chaine de caractere en entrée
		  perm_pos_init[(*perm_pos)][pos] = j;
		  perm_list[(*perm_pos)][pos] = seq_copy[j];
		  pos ++;
		}
	    (*perm_pos) ++; //Permutation computed !
	    //std::cout << (*perm_pos)  << std::endl;
	  }
	else if(n-nb_mask > k) //Pruning
	  {
	    get_perm_rec(seq_copy,n, perm_list,perm_pos_init,k,perm_pos,i+1,nb_mask+1,seq_size,null_element);
	  }
	delete [] seq_copy;
      }
  }

  static void swap_row(CImg<bool>& mat, int r1, int r2);
public:
  /*
   *Returns the number of permutations found
   * perm_list : permutations found ( allocation perfomed in function)
   * perm_pos_init : positions of the permutations in seq ( allocation perfomed in function)
   * s : T sequence to split
   * n   : sequence size
   * k : perm size
   */
  template<typename T>
  static int get_perm(const T * s, int nb_neighbours, T ** &perm_list, int ** &perm_pos_init, int k, T& null_element){
    int nb_perm =  gsl_sf_choose(nb_neighbours, k);
    /* Mem Alloc */
    perm_list = new T*[nb_perm];
    perm_pos_init = new int*[nb_perm];
    for(int i=0;i<nb_perm;i++)
      {
	perm_list[i] = new T[k];
	perm_pos_init[i] = new int[k];
      }
    //Cas particulier n = k
    if (nb_perm ==1)
      {
	for(int i=0;i<k;++i){
	  T tmp = s[i];
	  perm_list[0][i] = tmp;
	  //perm_list[0][i] = s[i];/*XXX: Invalid read of size sizeof(T)*/
	}
	for(int i=0;i<k;i++){
	  perm_pos_init[0][i] = i;
	}
      }
    else
      {
  	int perm_pos = 0;
  	T* seq = new T[nb_neighbours];
  	seq = (T*) memcpy(seq,s,nb_neighbours*sizeof(T));
  	get_perm_rec(seq,nb_neighbours,perm_list,perm_pos_init,k,&perm_pos,0,0,nb_neighbours, null_element);
	//delete []  seq; //[]
      }
    return nb_perm;
  }

  
  template<typename T>
  static int coucou(T i){
    return 45;
  }
  
  // template<typename T>
  // static int get_perm(const T * s, int nb_neighbours, T ** &perm_list, int ** &perm_pos_init, int k, T null_element){
    
  //   gsl_combination * c;
  //   size_t i;
  //   int nb_perm =  gsl_sf_choose(nb_neighbours, k);
  //   c = gsl_combination_calloc (nb_neighbours, k);

  //   perm_list = new T*[nb_perm];
  //   perm_pos_init = new int*[nb_perm];
  //   for(int i=0;i<nb_perm;i++)
  //     {
  // 	perm_pos_init[i] = new int[k];
  // 	perm_list[i] = new T[k];
  //     }
  //   int n = 0;
  //   do
  //     {
  // 	size_t * data = gsl_combination_data (c);
  // 	for(int i=0;i<k;i++)
  // 	  {
  // 	    perm_pos_init[n][i] = data[i];
  // 	    perm_list[n][i] = s[data[i]];
  // 	  }
  // 	n++;
  // 	//	gsl_combination_fprintf (stdout, c, " %u");
  //     }
  //   while (gsl_combination_next (c) == GSL_SUCCESS);
  //   gsl_combination_free (c);
    
  //   return nb_perm;
  // }
  static int booleanRank(CImg<bool> mat);
  
  static std::string reverse(std::string s);
  
  static int nbNeighbours(long i, pandore::Graph3d& g);


template<typename T>
 static void permute(T * tab,int pos1,int pos2) 
  { 
    T temp = tab[pos1];
    tab[pos1] = tab[pos2];
    tab[pos2] = temp;
  }

  /**
   * Sort list and return the number of permutation needed to do it.
   * The operator < must be define for T.
   */
 template<typename T>
 static  int countPermute(T * list,int size)
   {
     int nbPermutation=0;

     for(int i=0; i<size; i++)
       {
	 for(int j=i+1; j<size; j++)
	   {
	     if(list[i]<list[j])
	       {
		 nbPermutation++;
		 permute(list,i,j);
	       }
	   }
       }
     return nbPermutation;
   }

 /**
  * Return the position of e in tab l of size s (return -1 if e is not found)
  */
 template<typename T>
   static unsigned int pos(T * l,unsigned int s,T e)
   {
     for(unsigned int i=0;i<s;i++)
       if(l[i]==e) return i;
     return -1;
   }
 
};


#endif // __UTILS_H__

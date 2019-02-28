/*
 * @file utils.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Jun 13 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "utils.h"
#include <iostream>
using namespace std;

int utils::nbNeighbours(long i, pandore::Graph3d& g){
    int nb_neighbours = 0;
    for(pandore::GEdge * e = g[i]->Neighbours();e!=NULL;e=e->Next(),nb_neighbours++)
      ;
    return nb_neighbours;
}

string utils::reverse(string s)
{
  string reverse(s);
  for(unsigned int c=0;c<s.size();c++) 
    reverse[c] = s[s.size()-c-1];
  return reverse;
}



void utils::swap_row(CImg<bool>& mat, int r1, int r2)
{
  for(int i=0;i<mat.width();i++)
    {
      bool tmp = mat(i,r1);
      mat(i,r1) = mat(i,r2);
      mat(i,r2) = tmp;
    }  
}

int utils::booleanRank(CImg<bool> mat){
  int nr = mat.height();
  int nc = mat.width();
  
  for(int cur_row=0;cur_row<nr;cur_row++)
    {
      int cur_col;
      //Trouver le pivot
      for(int cur_row=0;cur_row<nr;cur_row++)
	{
	  for(int c=0;c<nc;c++)
	    {
	      for(int r=cur_row;r<nr;r++)
		if(mat(c,r))//Mon pivot
		  {
		    swap_row(mat,cur_row,r);	  
		    cur_col = c;
		  }
	    }
	  //Pivot = mat(cur_col, cur_row)
	  for(int r=cur_row+1;r<nr;r++)
	    for(int c=0;c<nc;c++)
	      mat(c,r) = mat(c,r)^mat(c,cur_row);
	}
    }
  
  int rank=0;
  for(int r=0;r<nr;r++)
    {
      int sum=0;
      for(int c=0;c<nc;c++)
	sum += mat(c,r);
      if(!sum)
	return rank;
      else
	rank++;
    }
  return rank;  
}

/******************************/
/*UTILS*/
/******************************/
// template<typename T>
// void utils::get_perm_rec(const T * seq, int n, T ** perm_list, int ** perm_pos_init, 
// 		  int k, int * perm_pos, int start_perm, int nb_mask, int seq_size, T null_element)
// {
//   for(int i=start_perm; i<n; i++)
//     {
//       T* seq_copy = new T[seq_size];
//       memcpy(seq_copy,seq,seq_size*sizeof(T));
//       seq_copy[i] = null_element;
//       if ((n-nb_mask-1) == k)     //Permutation #perm_pos get !
// 	{
// 	  int pos=0;
// 	  for(int j=0; j < seq_size; j++) 
// 	    if(seq_copy[j] != null_element)
// 	      {
// 		//j : position dans la chaine de caractere en entrée
// 		perm_pos_init[(*perm_pos)][pos] = j;
// 		perm_list[(*perm_pos)][pos] = seq_copy[j];
// 		pos ++;
// 	      }
// 	  (*perm_pos) ++; //Permutation computed !
// 	}
//       else if(n-nb_mask > k) //Pruning
// 	{
// 	  get_perm_rec(seq_copy,n, perm_list,perm_pos_init,k,perm_pos,i+1,nb_mask+1,seq_size,null_element);
// 	}
//       //XXX:
//       //delete [] seq_copy;
//     }
// }

// template<typename T>
// int utils::get_perm(const T * s, int nb_neighbours, 
// 	     T ** &perm_list, int ** &perm_pos_init, int k, T null_element)
// {
//   int nb_perm =  gsl_sf_choose(nb_neighbours, k);
//   /* Mem Alloc */
//   perm_list = new T*[nb_perm];
//   perm_pos_init = new int*[nb_perm];
//   for(int i=0;i<nb_perm;i++)
//     {
//       perm_pos_init[i] = new int[k];
//       perm_list[i] = new T[k];
//     }
//   //Cas particulier n = k
//   if (nb_perm ==1)
//     {
//       // (*perm_list)[0] = 
//       memcpy(perm_list[0],s,k*sizeof(T));
//       for(int i=0;i<k;i++)
//   	perm_pos_init[0][i] = i;
//     }
//   else
//     {
//       int perm_pos = 0;
//       T* seq = new T[nb_neighbours];
//       memcpy(seq,s,nb_neighbours*sizeof(T));
//       get_perm_rec(seq,nb_neighbours,perm_list,perm_pos_init,k,&perm_pos,0,0,nb_neighbours, null_element);
//     }
//   return nb_perm;
// }

// template<typename T> int utils::coucou(T i){
//   return 45;
// }


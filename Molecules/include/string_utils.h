/**
 * @file string_utils.h
 * @author Benoit Gauzere <<benoit.gauzere@ensicaen.fr>> 
 * @version     0.0.1 - Fri Feb  4 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __STRING_UTILS_H__
#define __STRING_UTILS_H__

#include <string>

class string_utils
{
private:
  static void get_perm_rec(std::string seq, int n, char ** perm_list, int ** perm_pos_init,
			   int k, int * perm_pos, int start_perm, int nb_mask);
public:
  static char * reverse_string(const char * orig);
  static bool keyComp(char * l, char * r);
  static bool keyStringComp(std::string l, std::string r);
  static char * stringSort(const char * to_sort);
  static int genericStringComp(const void * l, const void * r);
  /*
   *Returns the number of permutations found
   * perm_list : permutations found ( allocation perfomed in function)
   * perm_pos_init : positions of the permutations in seq ( allocation perfomed in function)
   * seq : char sequence to split
   * n   : sequence size
   * k : perm size
   */
  //TODO: A Templatiser <T = char>
  static int get_perm(char * seq, int seq_size, char *** perm_list,int *** perm_pos_init,  int k); 
};

#endif // __STRING_UTILS_H__

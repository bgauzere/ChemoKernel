/**
 * @file TreeletEditDistance.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon May 14 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __TREELETEDITDISTANCE_H__
#define __TREELETEDITDISTANCE_H__

#include "EditDistance.h"
#include <string>
#include <vector>



class TreeletEditDistance : public EditDistance
{
  static void *** edit_paths;
  static void *** edit_paths_edges;
  static void *** perm_intra_treelets;
  static int nb_perm_intra_treelets[];

  static int nb_perms[14][14]; //equals for edges
  static int size_perms[14][14]; //-1 for edges (a passer en [14]


  static const long G10toG0[6][1];
  static const long G10toG1[5][2];
  static const long G10toG2[5][3];
  static const long G10toG3[3][4];
  static const long G10toG4[2][5];
  static const long G10toG6[1][4];
  static const long G10toG7[1][5];
  static const long G11toG0[6][1];
  static const long G11toG1[5][2];
  static const long G11toG2[7][3];
  static const long G11toG3[3][4];
  static const long G11toG6[4][4];
  static const long G11toG7[3][5];
  static const long G11toG8[1][5];
  static const long G12toG0[6][1];
  static const long G12toG1[5][2];
  static const long G12toG2[6][3];
  static const long G12toG3[4][4];
  static const long G12toG6[2][4];
  static const long G12toG7[4][5];
  static const long G13toG0[6][1];
  static const long G13toG1[5][2];
  static const long G13toG2[10][3];
  static const long G13toG6[10][4];
  static const long G13toG8[5][5];
  static const long G1toG0[2][1];
  static const long G2toG0[3][1];
  static const long G2toG1[2][2];
  static const long G3toG0[4][1];
  static const long G3toG1[3][2];
  static const long G3toG2[2][3];
  static const long G4toG0[5][1];
  static const long G4toG1[4][2];
  static const long G4toG2[3][3];
  static const long G4toG3[2][4];
  static const long G5toG0[6][1];
  static const long G5toG1[5][2];
  static const long G5toG2[4][3];
  static const long G5toG3[3][4];
  static const long G5toG4[2][5];
  static const long G6toG0[4][1];
  static const long G6toG1[3][2];
  static const long G6toG2[3][3];
  static const long G7toG0[5][1];
  static const long G7toG1[4][2];
  static const long G7toG2[4][3];
  static const long G7toG3[2][4];
  static const long G7toG6[1][4];
  static const long G8toG0[5][1];
  static const long G8toG1[4][2];
  static const long G8toG2[6][3];
  static const long G8toG6[4][4];
  static const long G9toG0[6][1];
  static const long G9toG1[5][2];
  static const long G9toG2[5][3];
  static const long G9toG3[4][4];
  static const long G9toG4[1][5];
  static const long G9toG6[1][4];
  static const long G9toG7[2][5];
  static const long edgesG10toG1[5][1];
  static const long edgesG10toG2[5][2];
  static const long edgesG10toG3[3][3];
  static const long edgesG10toG4[2][4];
  static const long edgesG10toG6[1][3];
  static const long edgesG10toG7[1][4];
  static const long edgesG11toG1[5][1];
  static const long edgesG11toG2[7][2];
  static const long edgesG11toG3[3][3];
  static const long edgesG11toG6[4][3];
  static const long edgesG11toG7[3][4];
  static const long edgesG11toG8[1][4];
  static const long edgesG12toG1[5][1];
  static const long edgesG12toG2[6][2];
  static const long edgesG12toG3[4][3];
  static const long edgesG12toG6[2][3];
  static const long edgesG12toG7[4][4];
  static const long edgesG13toG1[5][1];
  static const long edgesG13toG2[10][2];
  static const long edgesG13toG6[10][3];
  static const long edgesG13toG8[5][4];
  static const long edgesG2toG1[2][1];
  static const long edgesG3toG1[3][1];
  static const long edgesG3toG2[2][2];
  static const long edgesG4toG1[4][1];
  static const long edgesG4toG2[3][2];
  static const long edgesG4toG3[2][3];
  static const long edgesG5toG1[5][1];
  static const long edgesG5toG2[4][2];
  static const long edgesG5toG3[3][3];
  static const long edgesG5toG4[2][4];
  static const long edgesG6toG1[3][1];
  static const long edgesG6toG2[3][2];
  static const long edgesG7toG1[4][1];
  static const long edgesG7toG2[4][2];
  static const long edgesG7toG3[2][3];
  static const long edgesG7toG6[1][3];
  static const long edgesG8toG1[4][1];
  static const long edgesG8toG2[6][2];
  static const long edgesG8toG6[4][3];
  static const long edgesG9toG1[5][1];
  static const long edgesG9toG2[5][2];
  static const long edgesG9toG3[4][3];
  static const long edgesG9toG4[1][4];
  static const long edgesG9toG6[1][3];
  static const long edgesG9toG7[2][4];
  

  /*TODO : Permutations Intra-Treelets
    static const int permutations_treelets[14][]*/
  
  //#include "TreeletEditDistanceProjections.h" 
  static double edit_distance_struct[14][14];
  static double mcs[14][14];
  static double StructuralCost(int treelet_i, int treelet_j);
  static long offsetNode(int perm,int pos,int ti,int tj);
  static long offsetEdge(int perm,int pos,int ti,int tj);
  double substitution_cost;
  double deletion_cost;
public:
  TreeletEditDistance(double c_s,double c_d);
  static void Init();
  static std::vector<std::string> Codes(int treelet, std::string code,int treelet_dest);
  static int compareNodesCodes(std::string code_i,std::string code_j,int treelet);
  static int maxPermsToCompute();
  double operator() (int t_i, std::string c_i, int t_j, std::string c_j);
  double operator() ( pandore::Collection * c1, pandore::Collection * c2 );
};


#endif // __TREELETEDITDISTANCE_H__



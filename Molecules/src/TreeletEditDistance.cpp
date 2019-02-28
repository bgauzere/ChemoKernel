/*
 * @file TreeletEditDistance.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon May 14 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "TreeletEditDistance.h"
#include <iostream>
#include <string.h>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include <limits.h>
#include <float.h>
using namespace std;
using namespace cimg_library;

void *** TreeletEditDistance::edit_paths=NULL;
void *** TreeletEditDistance::edit_paths_edges=NULL;
void *** TreeletEditDistance::perm_intra_treelets=NULL;

double TreeletEditDistance::edit_distance_struct[14][14] = {
  {0,1,2,3,4,5,3,4,4,5,5,5,5,5},
  {1,0,1,2,3,4,2,3,3,4,4,4,4,4},
  {2,1,0,1,2,3,1,2,2,3,3,3,3,3},
  {3,2,1,0,1,2,2,1,3,2,2,2,2,4},
  {4,3,2,1,0,1,3,2,4,1,1,3,3,5},
  {5,4,3,2,1,0,4,3,5,2,2,4,4,6},
  {3,2,1,2,3,4,0,1,1,2,2,2,2,2},
  {4,3,2,1,2,3,1,0,2,1,1,1,1,3},
  {4,3,2,3,4,5,1,2,0,3,3,1,3,1},
  {5,4,3,2,1,2,2,1,3,0,2,2,2,4},
  {5,4,3,2,1,2,2,1,3,2,0,2,2,4},
  {5,4,3,2,3,4,2,1,1,2,2,0,2,2},
  {5,4,3,2,3,4,2,1,3,2,2,2,0,4},
  {5,4,3,4,5,6,2,3,1,4,4,2,4,0}
};
double TreeletEditDistance::mcs[14][14] = {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {0,1,1,1,1,1,1,1,1,1,1,1,1,1},
  {0,1,2,2,2,2,2,2,2,2,2,2,2,2},
  {0,1,2,3,3,3,2,3,2,3,3,3,3,2},
  {0,1,2,3,4,4,2,3,2,4,4,3,3,2},
  {0,1,2,3,4,5,2,3,2,4,4,3,3,2},
  {0,1,2,2,2,2,6,6,6,6,6,6,6,6},
  {0,1,2,3,3,3,6,7,6,7,7,7,7,6},
  {0,1,2,2,2,2,6,6,8,6,6,8,6,8},
  {0,1,2,3,4,4,6,7,6,9,7,7,7,6},
  {0,1,2,3,4,4,6,7,6,7,10,7,7,6},
  {0,1,2,3,3,3,6,7,8,7,7,11,7,8},
  {0,1,2,3,3,3,6,7,6,7,7,7,12,6},
  {0,1,2,2,2,2,6,6,8,6,6,8,6,13}
};

int TreeletEditDistance::nb_perms[14][14]={
  {1,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {2,1,0,0,0,0,0,0,0,0,0,0,0,0},
  {3,2,1,0,0,0,0,0,0,0,0,0,0,0},
  {4,3,2,1,0,0,0,0,0,0,0,0,0,0},
  {5,4,3,2,1,0,0,0,0,0,0,0,0,0},
  {6,5,4,3,2,1,0,0,0,0,0,0,0,0},
  {4,3,3,0,0,0,1,0,0,0,0,0,0,0},
  {5,4,4,2,0,0,1,1,0,0,0,0,0,0},
  {5,4,6,0,0,0,4,0,1,0,0,0,0,0},
  {6,5,5,4,1,0,1,2,0,1,0,0,0,0},
  {6,5,5,3,2,0,1,1,0,0,1,0,0,0},
  {6,5,7,3,0,0,4,3,1,0,0,1,0,0},
  {6,5,6,4,0,0,2,4,0,0,0,0,1,0},
  {6,5,10,0,0,0,10,0,5,0,0,0,0,1}
};

int TreeletEditDistance::size_perms[14][14]={
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {1,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {1,2,0,0,0,0,0,0,0,0,0,0,0,0},
  {1,2,3,0,0,0,0,0,0,0,0,0,0,0},
  {1,2,3,4,0,0,0,0,0,0,0,0,0,0},
  {1,2,3,4,5,0,0,0,0,0,0,0,0,0},
  {1,2,3,0,0,0,0,0,0,0,0,0,0,0},
  {1,2,3,4,0,0,4,0,0,0,0,0,0,0},
  {1,2,3,0,0,0,4,0,0,0,0,0,0,0},
  {1,2,3,4,5,0,4,5,0,0,0,0,0,0},
  {1,2,3,4,5,0,4,5,0,0,0,0,0,0},
  {1,2,3,4,0,0,4,5,5,0,0,0,0,0},
  {1,2,3,4,0,0,4,5,0,0,0,0,0,0},
  {1,2,3,0,0,0,4,0,5,0,0,0,0,0}
};

int TreeletEditDistance::maxPermsToCompute(){
  int max =0;
  for(int i=0;i<14;i++)
    for(int j=i;j<14;j++)
      {//Calcul de n_i*n_j
	int mcs_test =  TreeletEditDistance::mcs[i][j];
       	int nb_test = nb_perms[i][mcs_test]*nb_perms[j][mcs_test];
	if (i ==j) 
	  nb_test =1;
	
	if ((mcs_test<6) && (mcs_test > 1))
	  nb_test=nb_test*2;
	if ((mcs_test==6)||(mcs_test==11))
	  nb_test=nb_test*6;
	if((mcs_test==7) || (mcs_test==9) || (mcs_test==10))
	  nb_test=nb_test*2;
	if(mcs_test==8)
	  nb_test=nb_test*24;
	if(mcs_test==12)
	  nb_test=nb_test*8;
	if(mcs_test==13)
	  nb_test=nb_test*120;
	cout << "(" << i << ","<< j << ")," << mcs_test << ","<< nb_test << endl;
	max=(nb_test>max)?nb_test:max;
      }
  return max;
}
TreeletEditDistance::TreeletEditDistance(double c_s,double c_d){
  substitution_cost = c_s;
  deletion_cost = c_d;  
}

double TreeletEditDistance::StructuralCost(int i, int j){
  return edit_distance_struct[i][j];
}


void TreeletEditDistance::Init(){

  /*Iniatialisation des permutations*/
  TreeletEditDistance::perm_intra_treelets = new void**[14];
  for(int i=0;i<14;i++){
    TreeletEditDistance::perm_intra_treelets[i] = new void*[14];
    memset(perm_intra_treelets[i],0,sizeof(void*)*14);
  }
  
  TreeletEditDistance::edit_paths = new void**[14];
  for(int i=0;i<14;i++){
    TreeletEditDistance::edit_paths[i] = new void*[14];
    memset(edit_paths[i],0,sizeof(void*)*14);
  }
  

  TreeletEditDistance::edit_paths[1][0] = (void*)(TreeletEditDistance::G1toG0);
  TreeletEditDistance::edit_paths[2][0] = (void*)(TreeletEditDistance::G2toG0);
  TreeletEditDistance::edit_paths[2][1] = (void*)(TreeletEditDistance::G2toG1);
  TreeletEditDistance::edit_paths[3][0] = (void*)(TreeletEditDistance::G3toG0);
  TreeletEditDistance::edit_paths[3][1] = (void*)(TreeletEditDistance::G3toG1);
  TreeletEditDistance::edit_paths[3][2] = (void*)(TreeletEditDistance::G3toG2);
  TreeletEditDistance::edit_paths[4][0] = (void*)(TreeletEditDistance::G4toG0);
  TreeletEditDistance::edit_paths[4][1] = (void*)(TreeletEditDistance::G4toG1);
  TreeletEditDistance::edit_paths[4][2] = (void*)(TreeletEditDistance::G4toG2);
  TreeletEditDistance::edit_paths[4][3] = (void*)(TreeletEditDistance::G4toG3);
  TreeletEditDistance::edit_paths[5][0] = (void*)(TreeletEditDistance::G5toG0);
  TreeletEditDistance::edit_paths[5][1] = (void*)(TreeletEditDistance::G5toG1);
  TreeletEditDistance::edit_paths[5][2] = (void*)(TreeletEditDistance::G5toG2);
  TreeletEditDistance::edit_paths[5][3] = (void*)(TreeletEditDistance::G5toG3);
  TreeletEditDistance::edit_paths[5][4] = (void*)(TreeletEditDistance::G5toG4);
  TreeletEditDistance::edit_paths[6][0] = (void*)(TreeletEditDistance::G6toG0);
  TreeletEditDistance::edit_paths[6][1] = (void*)(TreeletEditDistance::G6toG1);
  TreeletEditDistance::edit_paths[6][2] = (void*)(TreeletEditDistance::G6toG2);
  TreeletEditDistance::edit_paths[7][0] = (void*)(TreeletEditDistance::G7toG0);
  TreeletEditDistance::edit_paths[7][1] = (void*)(TreeletEditDistance::G7toG1);
  TreeletEditDistance::edit_paths[7][2] = (void*)(TreeletEditDistance::G7toG2);
  TreeletEditDistance::edit_paths[7][3] = (void*)(TreeletEditDistance::G7toG3);
  TreeletEditDistance::edit_paths[7][6] = (void*)(TreeletEditDistance::G7toG6);
  TreeletEditDistance::edit_paths[8][0] = (void*)(TreeletEditDistance::G8toG0);
  TreeletEditDistance::edit_paths[8][1] = (void*)(TreeletEditDistance::G8toG1);
  TreeletEditDistance::edit_paths[8][2] = (void*)(TreeletEditDistance::G8toG2);
  TreeletEditDistance::edit_paths[8][6] = (void*)(TreeletEditDistance::G8toG6);
  TreeletEditDistance::edit_paths[9][0] = (void*)(TreeletEditDistance::G9toG0);
  TreeletEditDistance::edit_paths[9][1] = (void*)(TreeletEditDistance::G9toG1);
  TreeletEditDistance::edit_paths[9][2] = (void*)(TreeletEditDistance::G9toG2);
  TreeletEditDistance::edit_paths[9][3] = (void*)(TreeletEditDistance::G9toG3);
  TreeletEditDistance::edit_paths[9][4] = (void*)(TreeletEditDistance::G9toG4);
  TreeletEditDistance::edit_paths[9][6] = (void*)(TreeletEditDistance::G9toG6);
  TreeletEditDistance::edit_paths[9][7] = (void*)(TreeletEditDistance::G9toG7);
  TreeletEditDistance::edit_paths[10][0] = (void*)(TreeletEditDistance::G10toG0);
  TreeletEditDistance::edit_paths[10][1] = (void*)(TreeletEditDistance::G10toG1);
  TreeletEditDistance::edit_paths[10][2] = (void*)(TreeletEditDistance::G10toG2);
  TreeletEditDistance::edit_paths[10][3] = (void*)(TreeletEditDistance::G10toG3);
  TreeletEditDistance::edit_paths[10][4] = (void*)(TreeletEditDistance::G10toG4);
  TreeletEditDistance::edit_paths[10][6] = (void*)(TreeletEditDistance::G10toG6);
  TreeletEditDistance::edit_paths[10][7] = (void*)(TreeletEditDistance::G10toG7);
  TreeletEditDistance::edit_paths[11][0] = (void*)(TreeletEditDistance::G11toG0);
  TreeletEditDistance::edit_paths[11][1] = (void*)(TreeletEditDistance::G11toG1);
  TreeletEditDistance::edit_paths[11][2] = (void*)(TreeletEditDistance::G11toG2);
  TreeletEditDistance::edit_paths[11][3] = (void*)(TreeletEditDistance::G11toG3);
  TreeletEditDistance::edit_paths[11][6] = (void*)(TreeletEditDistance::G11toG6);
  TreeletEditDistance::edit_paths[11][7] = (void*)(TreeletEditDistance::G11toG7);
  TreeletEditDistance::edit_paths[11][8] = (void*)(TreeletEditDistance::G11toG8);
  TreeletEditDistance::edit_paths[12][0] = (void*)(TreeletEditDistance::G12toG0);
  TreeletEditDistance::edit_paths[12][1] = (void*)(TreeletEditDistance::G12toG1);
  TreeletEditDistance::edit_paths[12][2] = (void*)(TreeletEditDistance::G12toG2);
  TreeletEditDistance::edit_paths[12][3] = (void*)(TreeletEditDistance::G12toG3);
  TreeletEditDistance::edit_paths[12][6] = (void*)(TreeletEditDistance::G12toG6);
  TreeletEditDistance::edit_paths[12][7] = (void*)(TreeletEditDistance::G12toG7);
  TreeletEditDistance::edit_paths[13][0] = (void*)(TreeletEditDistance::G13toG0);
  TreeletEditDistance::edit_paths[13][1] = (void*)(TreeletEditDistance::G13toG1);
  TreeletEditDistance::edit_paths[13][2] = (void*)(TreeletEditDistance::G13toG2);
  TreeletEditDistance::edit_paths[13][6] = (void*)(TreeletEditDistance::G13toG6);
  TreeletEditDistance::edit_paths[13][8] = (void*)(TreeletEditDistance::G13toG8);

  TreeletEditDistance::edit_paths_edges = new void**[14];
  for(int i=0;i<14;i++){
    TreeletEditDistance::edit_paths_edges[i] = new void*[14];
    memset(edit_paths_edges[i],0,sizeof(void*)*14);
  }
  
  TreeletEditDistance::edit_paths_edges[2][1] = (void*)(TreeletEditDistance::edgesG2toG1);
  TreeletEditDistance::edit_paths_edges[3][1] = (void*)(TreeletEditDistance::edgesG3toG1);
  TreeletEditDistance::edit_paths_edges[3][2] = (void*)(TreeletEditDistance::edgesG3toG2);
  TreeletEditDistance::edit_paths_edges[4][1] = (void*)(TreeletEditDistance::edgesG4toG1);
  TreeletEditDistance::edit_paths_edges[4][2] = (void*)(TreeletEditDistance::edgesG4toG2);
  TreeletEditDistance::edit_paths_edges[4][3] = (void*)(TreeletEditDistance::edgesG4toG3);
  TreeletEditDistance::edit_paths_edges[5][1] = (void*)(TreeletEditDistance::edgesG5toG1);
  TreeletEditDistance::edit_paths_edges[5][2] = (void*)(TreeletEditDistance::edgesG5toG2);
  TreeletEditDistance::edit_paths_edges[5][3] = (void*)(TreeletEditDistance::edgesG5toG3);
  TreeletEditDistance::edit_paths_edges[5][4] = (void*)(TreeletEditDistance::edgesG5toG4);
  TreeletEditDistance::edit_paths_edges[6][1] = (void*)(TreeletEditDistance::edgesG6toG1);
  TreeletEditDistance::edit_paths_edges[6][2] = (void*)(TreeletEditDistance::edgesG6toG2);
  TreeletEditDistance::edit_paths_edges[7][1] = (void*)(TreeletEditDistance::edgesG7toG1);
  TreeletEditDistance::edit_paths_edges[7][2] = (void*)(TreeletEditDistance::edgesG7toG2);
  TreeletEditDistance::edit_paths_edges[7][3] = (void*)(TreeletEditDistance::edgesG7toG3);
  TreeletEditDistance::edit_paths_edges[7][6] = (void*)(TreeletEditDistance::edgesG7toG6);
  TreeletEditDistance::edit_paths_edges[8][1] = (void*)(TreeletEditDistance::edgesG8toG1);
  TreeletEditDistance::edit_paths_edges[8][2] = (void*)(TreeletEditDistance::edgesG8toG2);
  TreeletEditDistance::edit_paths_edges[8][6] = (void*)(TreeletEditDistance::edgesG8toG6);
  TreeletEditDistance::edit_paths_edges[9][1] = (void*)(TreeletEditDistance::edgesG9toG1);
  TreeletEditDistance::edit_paths_edges[9][2] = (void*)(TreeletEditDistance::edgesG9toG2);
  TreeletEditDistance::edit_paths_edges[9][3] = (void*)(TreeletEditDistance::edgesG9toG3);
  TreeletEditDistance::edit_paths_edges[9][4] = (void*)(TreeletEditDistance::edgesG9toG4);
  TreeletEditDistance::edit_paths_edges[9][6] = (void*)(TreeletEditDistance::edgesG9toG6);
  TreeletEditDistance::edit_paths_edges[9][7] = (void*)(TreeletEditDistance::edgesG9toG7);
  TreeletEditDistance::edit_paths_edges[10][1] = (void*)(TreeletEditDistance::edgesG10toG1);
  TreeletEditDistance::edit_paths_edges[10][2] = (void*)(TreeletEditDistance::edgesG10toG2);
  TreeletEditDistance::edit_paths_edges[10][3] = (void*)(TreeletEditDistance::edgesG10toG3);
  TreeletEditDistance::edit_paths_edges[10][4] = (void*)(TreeletEditDistance::edgesG10toG4);
  TreeletEditDistance::edit_paths_edges[10][6] = (void*)(TreeletEditDistance::edgesG10toG6);
  TreeletEditDistance::edit_paths_edges[10][7] = (void*)(TreeletEditDistance::edgesG10toG7);
  TreeletEditDistance::edit_paths_edges[11][1] = (void*)(TreeletEditDistance::edgesG11toG1);
  TreeletEditDistance::edit_paths_edges[11][2] = (void*)(TreeletEditDistance::edgesG11toG2);
  TreeletEditDistance::edit_paths_edges[11][3] = (void*)(TreeletEditDistance::edgesG11toG3);
  TreeletEditDistance::edit_paths_edges[11][6] = (void*)(TreeletEditDistance::edgesG11toG6);
  TreeletEditDistance::edit_paths_edges[11][7] = (void*)(TreeletEditDistance::edgesG11toG7);
  TreeletEditDistance::edit_paths_edges[11][8] = (void*)(TreeletEditDistance::edgesG11toG8);
  TreeletEditDistance::edit_paths_edges[12][1] = (void*)(TreeletEditDistance::edgesG12toG1);
  TreeletEditDistance::edit_paths_edges[12][2] = (void*)(TreeletEditDistance::edgesG12toG2);
  TreeletEditDistance::edit_paths_edges[12][3] = (void*)(TreeletEditDistance::edgesG12toG3);
  TreeletEditDistance::edit_paths_edges[12][6] = (void*)(TreeletEditDistance::edgesG12toG6);
  TreeletEditDistance::edit_paths_edges[12][7] = (void*)(TreeletEditDistance::edgesG12toG7);
  TreeletEditDistance::edit_paths_edges[13][1] = (void*)(TreeletEditDistance::edgesG13toG1);
  TreeletEditDistance::edit_paths_edges[13][2] = (void*)(TreeletEditDistance::edgesG13toG2);
  TreeletEditDistance::edit_paths_edges[13][6] = (void*)(TreeletEditDistance::edgesG13toG6);
  TreeletEditDistance::edit_paths_edges[13][8] = (void*)(TreeletEditDistance::edgesG13toG8);
  
}

/*Node Permutations*/
const long TreeletEditDistance::G1toG0[2][1] = {{0},{1}};
const long TreeletEditDistance::G2toG0[3][1] = {{0},{1},{2}};
const long TreeletEditDistance::G2toG1[2][2] = {{0,1},{1,2}};
const long TreeletEditDistance::G3toG0[4][1] = {{0},{1},{2},{3}};
const long TreeletEditDistance::G3toG1[3][2] = {{0,1},{1,2},{2,3}};
const long TreeletEditDistance::G3toG2[2][3] = {{0,1,2},{1,2,3}};
const long TreeletEditDistance::G4toG0[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::G4toG1[4][2] = {{0,1},{1,2},{2,3},{3,4}};
const long TreeletEditDistance::G4toG2[3][3] = {{0,1,2},{1,2,3},{2,3,4}};
const long TreeletEditDistance::G4toG3[2][4] = {{0,1,2,3},{1,2,3,4}};
const long TreeletEditDistance::G5toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G5toG1[5][2] = {{0,1},{1,2},{2,3},{3,4},{4,5}};
const long TreeletEditDistance::G5toG2[4][3] = {{0,1,2},{1,2,3},{2,3,4},{3,4,5}};
const long TreeletEditDistance::G5toG3[3][4] = {{0,1,2,3},{1,2,3,4},{2,3,4,5}};
const long TreeletEditDistance::G5toG4[2][5] = {{0,1,2,3,4},{1,2,3,4,5}};
const long TreeletEditDistance::G6toG0[4][1] = {{0},{1},{2},{3}};
const long TreeletEditDistance::G6toG1[3][2] = {{0,1},{0,2},{0,3}};
const long TreeletEditDistance::G6toG2[3][3] = {{1,0,2},{1,0,3},{2,0,3}};
const long TreeletEditDistance::G7toG0[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::G7toG1[4][2] = {{0,1},{1,2},{2,3},{2,4}};
const long TreeletEditDistance::G7toG2[4][3] = {{0,1,2},{1,2,3},{1,2,4},{3,2,4}};
const long TreeletEditDistance::G7toG3[2][4] = {{0,1,2,3},{0,1,2,4}};
const long TreeletEditDistance::G7toG6[1][4] = {{2,1,3,4}};
const long TreeletEditDistance::G8toG0[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::G8toG1[4][2] = {{0,1},{0,2},{0,3},{0,4}};
const long TreeletEditDistance::G8toG2[6][3] = {{1,0,2},{1,0,3},{1,0,4},{2,0,3},{2,0,4},{3,0,4}};
const long TreeletEditDistance::G8toG6[4][4] = {{0,1,2,3},{0,1,2,4},{0,1,3,4},{0,2,3,4}};
const long TreeletEditDistance::G9toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G9toG1[5][2] = {{0,1},{0,2},{0,4},{2,3},{4,5}};
const long TreeletEditDistance::G9toG2[5][3] = {{0,2,3},{0,4,5},{1,0,2},{1,0,4},{2,0,4}};
const long TreeletEditDistance::G9toG3[4][4] = {{1,0,2,3},{1,0,4,5},{2,0,4,5},{3,2,0,4}};
const long TreeletEditDistance::G9toG4[1][5] = {{3,2,0,4,5}};
const long TreeletEditDistance::G9toG6[1][4] = {{0,1,2,4}};
const long TreeletEditDistance::G9toG7[2][5] = {{3,2,0,1,4},{5,4,0,1,2}};
const long TreeletEditDistance::G10toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G10toG1[5][2] = {{0,1},{1,2},{2,3},{3,4},{3,5}};
const long TreeletEditDistance::G10toG2[5][3] = {{0,1,2},{1,2,3},{2,3,4},{2,3,5},{4,3,5}};
const long TreeletEditDistance::G10toG3[3][4] = {{0,1,2,3},{1,2,3,4},{1,2,3,5}};
const long TreeletEditDistance::G10toG4[2][5] = {{0,1,2,3,4},{0,1,2,3,5}};
const long TreeletEditDistance::G10toG6[1][4] = {{3,2,4,5}};
const long TreeletEditDistance::G10toG7[1][5] = {{1,2,3,4,5}};
const long TreeletEditDistance::G11toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G11toG1[5][2] = {{0,1},{1,2},{2,3},{2,4},{2,5}};
const long TreeletEditDistance::G11toG2[7][3] = {{0,1,2},{1,2,3},{1,2,4},{1,2,5},{3,2,4},{3,2,5},{4,2,5}};
const long TreeletEditDistance::G11toG3[3][4] = {{0,1,2,3},{0,1,2,4},{0,1,2,5}};
const long TreeletEditDistance::G11toG6[4][4] = {{2,1,3,4},{2,1,3,5},{2,1,4,5},{2,3,4,5}};
const long TreeletEditDistance::G11toG7[3][5] = {{0,1,2,3,4},{0,1,2,3,5},{0,1,2,4,5}};
const long TreeletEditDistance::G11toG8[1][5] = {{2,1,3,4,5}};
const long TreeletEditDistance::G12toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G12toG1[5][2] = {{0,1},{0,2},{0,3},{3,4},{3,5}};
const long TreeletEditDistance::G12toG2[6][3] = {{0,3,4},{0,3,5},{1,0,2},{1,0,3},{2,0,3},{4,3,5}};
const long TreeletEditDistance::G12toG3[4][4] = {{1,0,3,4},{1,0,3,5},{2,0,3,4},{2,0,3,5}};
const long TreeletEditDistance::G12toG6[2][4] = {{0,1,2,3},{3,0,4,5}};
const long TreeletEditDistance::G12toG7[4][5] = {{1,0,3,4,5},{2,0,3,4,5},{4,3,0,1,2},{5,3,0,1,2}};
const long TreeletEditDistance::G13toG0[6][1] = {{0},{1},{2},{3},{4},{5}};
const long TreeletEditDistance::G13toG1[5][2] = {{0,1},{0,2},{0,3},{0,4},{0,5}};
const long TreeletEditDistance::G13toG2[10][3] = {{1,0,2},{1,0,3},{1,0,4},{1,0,5},{2,0,3},{2,0,4},{2,0,5},{3,0,4},{3,0,5},{4,0,5}};
const long TreeletEditDistance::G13toG6[10][4] = {{0,1,2,3},{0,1,2,4},{0,1,2,5},{0,1,3,4},{0,1,3,5},{0,1,4,5},{0,2,3,4},{0,2,3,5},{0,2,4,5},{0,3,4,5}};
const long TreeletEditDistance::G13toG8[5][5] = {{0,1,2,3,4},{0,1,2,3,5},{0,1,2,4,5},{0,1,3,4,5},{0,2,3,4,5}};

/*Edges Permutations*/
const long TreeletEditDistance::edgesG2toG1[2][1] = {{1},{1}};
const long TreeletEditDistance::edgesG3toG1[3][1] = {{0},{1},{2}};
const long TreeletEditDistance::edgesG3toG2[2][2] = {{0,1},{1,2}};
const long TreeletEditDistance::edgesG4toG1[4][1] = {{0},{1},{2},{3}};
const long TreeletEditDistance::edgesG4toG2[3][2] = {{0,1},{1,2},{2,3}};
const long TreeletEditDistance::edgesG4toG3[2][3] = {{0,1,2},{1,2,3}};
const long TreeletEditDistance::edgesG5toG1[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::edgesG5toG2[4][2] = {{0,1},{1,2},{2,3},{3,4}};
const long TreeletEditDistance::edgesG5toG3[3][3] = {{0,1,2},{1,2,3},{2,3,4}};
const long TreeletEditDistance::edgesG5toG4[2][4] = {{0,1,2,3},{1,2,3,4}};
const long TreeletEditDistance::edgesG6toG1[3][1] = {{0},{1},{2}};
const long TreeletEditDistance::edgesG6toG2[3][2] = {{0,1},{0,2},{1,2}};
const long TreeletEditDistance::edgesG7toG1[4][1] = {{0},{1},{2},{3}};
const long TreeletEditDistance::edgesG7toG2[4][2] = {{0,1},{1,2},{1,3},{2,3}};
const long TreeletEditDistance::edgesG7toG3[2][3] = {{0,1,2},{0,1,3}};
const long TreeletEditDistance::edgesG7toG6[1][3] = {{1,2,3}};
const long TreeletEditDistance::edgesG8toG1[4][1] = {{0},{1},{2},{3}};
const long TreeletEditDistance::edgesG8toG2[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
const long TreeletEditDistance::edgesG8toG6[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
const long TreeletEditDistance::edgesG9toG1[5][1] = {{0},{1},{3},{2},{4}};
const long TreeletEditDistance::edgesG9toG2[5][2] = {{1,2},{3,4},{0,1},{0,3},{1,3}};
const long TreeletEditDistance::edgesG9toG3[4][3] = {{0,1,2},{0,3,4},{1,3,4},{2,1,3}};
const long TreeletEditDistance::edgesG9toG4[1][4] = {{2,1,3,4}};
const long TreeletEditDistance::edgesG9toG6[1][3] = {{0,1,3}};
const long TreeletEditDistance::edgesG9toG7[2][4] = {{2,1,0,3},{4,3,0,1}};
const long TreeletEditDistance::edgesG10toG1[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::edgesG10toG2[5][2] = {{0,1},{1,2},{2,3},{2,4},{3,4}};
const long TreeletEditDistance::edgesG10toG3[3][3] = {{0,1,2},{1,2,3},{1,2,4}};
const long TreeletEditDistance::edgesG10toG4[2][4] = {{0,1,2,3},{0,1,2,4}};
const long TreeletEditDistance::edgesG10toG6[1][3] = {{2,3,4}};
const long TreeletEditDistance::edgesG10toG7[1][4] = {{1,2,3,4}};
const long TreeletEditDistance::edgesG11toG1[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::edgesG11toG2[7][2] = {{0,1},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
const long TreeletEditDistance::edgesG11toG3[3][3] = {{0,1,2},{0,1,3},{0,1,4}};
const long TreeletEditDistance::edgesG11toG6[4][3] = {{1,2,3},{1,2,4},{1,3,4},{2,3,4}};
const long TreeletEditDistance::edgesG11toG7[3][4] = {{0,1,2,3},{0,1,2,4},{0,1,3,4}};
const long TreeletEditDistance::edgesG11toG8[1][4] = {{1,2,3,4}};
const long TreeletEditDistance::edgesG12toG1[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::edgesG12toG2[6][2] = {{2,3},{2,4},{0,1},{0,2},{1,2},{3,4}};
const long TreeletEditDistance::edgesG12toG3[4][3] = {{0,2,3},{0,2,4},{1,2,3},{1,2,4}};
const long TreeletEditDistance::edgesG12toG6[2][3] = {{0,1,2},{2,3,4}};
const long TreeletEditDistance::edgesG12toG7[4][4] = {{0,2,3,4},{1,2,3,4},{3,2,0,1},{4,2,0,1}};
const long TreeletEditDistance::edgesG13toG1[5][1] = {{0},{1},{2},{3},{4}};
const long TreeletEditDistance::edgesG13toG2[10][2] = {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
const long TreeletEditDistance::edgesG13toG6[10][3] = {{0,1,2},{0,1,3},{0,1,4},{0,2,3},{0,2,4},{0,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4}};
const long TreeletEditDistance::edgesG13toG8[5][4] = {{0,1,2,3},{0,1,2,4},{0,1,3,4},{0,2,3,4},{1,2,3,4}};

long TreeletEditDistance::offsetNode(int perm,int pos,int ti,int tj){
  return(((long*)(TreeletEditDistance::edit_paths[ti][tj]))[perm*size_perms[ti][tj]+pos]);
}
long TreeletEditDistance::offsetEdge(int perm,int pos,int ti,int tj){
  return(((long*)(TreeletEditDistance::edit_paths_edges[ti][tj]))[perm*(size_perms[ti][tj]-1)+pos]);
}

vector<string> TreeletEditDistance::Codes(int treelet, std::string code,int treelet_dest){
  if(treelet == treelet_dest)
    return vector<string>(1,code);
  /*Node/edge separated by '_'*/
  vector<string> nodes;
  vector<string> edges;
  bool is_node = true;
  for(int c=0;c<code.length();c++)
    if(code[c]!='_'){
      if(is_node)
	nodes.push_back(string(1,code[c]));
      else
	edges.push_back(string(1,code[c]));
      is_node =! is_node;
    }
  assert (nb_perms[treelet][treelet_dest] != 0);  
  long ** tmp = (long**) edit_paths[treelet][treelet_dest];
  vector<string> codes_perms;
  for(int i=0;i<nb_perms[treelet][treelet_dest];i++){
    stringstream ss;
    for(int j=0;j<size_perms[treelet][treelet_dest];j++){
      ss << nodes[offsetNode(i,j,treelet,treelet_dest)];
      if(j!=size_perms[treelet][treelet_dest]-1)//#edges = #nodes -1
	ss<<"_"<<edges[offsetEdge(i,j,treelet,treelet_dest)]<<"_";
    }
    codes_perms.push_back(ss.str());
  }
  return codes_perms;
}

int  TreeletEditDistance::compareNodesCodes(string code_i,string code_j, int treelet_type){
  if (treelet_type < 6){//Linear Treelet
    int diff =0;
    for(int i=0;i<code_i.length();i++)
      diff += code_i[i] !=  code_j[i];
    int diff_2 =0;
    for(int i=0;i<code_i.length();i++)
      diff_2 += code_i[i] !=  code_j[code_j.length()-i-1];
    return (diff<diff_2)?diff:diff_2;
  }
  else if ((treelet_type == 6) || (treelet_type == 8) || (treelet_type == 13)){ //n-star
    int diff = (code_i[0] != code_j[0]);
    int nb_branches;
    switch(treelet_type){
    case 6:
      nb_branches = 3;
      break;
    case 8:
      nb_branches = 4;
      break;
    case 13:
      nb_branches = 5;
      break;
    }
    vector<pair<char,char> > branches_i(nb_branches);
    vector<pair<char,char> > branches_j(nb_branches);
    for(int i=0;i<nb_branches;i++){
      branches_i[i] = pair<char,char>(code_i[2+(i*4)],code_i[2+(i*4)+2]);
      branches_j[i] = pair<char,char>(code_j[2+(i*4)],code_j[2+(i*4)+2]);
    }
    CImg<double> cost(nb_branches,nb_branches);
    for(int i=0;i<nb_branches;i++)
      for(int j=0;j<nb_branches;j++)
	cost(i,j) = ((branches_i[i].first != branches_j[j].first) +
		     (branches_i[i].second != branches_j[j].second));
    Alignement * align = EditDistance::applyMunkres(cost);
    double val =diff + align->cost;
    delete align;
    return val;
    
  }else if((treelet_type == 7) || (treelet_type == 10)){//permut between e3-v4 and e4-v5
    int diff = 0;
    for(int i = 0; i<code_i.length()-(4*2);i++)
      diff += code_i[i] != code_j[i];
    vector<pair<char,char> > branches_i(2);
    vector<pair<char,char> > branches_j(2);
    branches_i[0] = pair<char,char>(code_i[code_i.length()-7],code_i[code_i.length()-5]);
    branches_i[1] = pair<char,char>(code_i[code_i.length()-3],code_i[code_i.length()-1]);
    branches_j[0] = pair<char,char>(code_j[code_j.length()-7],code_j[code_j.length()-5]);
    branches_j[1] = pair<char,char>(code_j[code_j.length()-3],code_j[code_j.length()-1]);
    CImg<double> cost(2,2);
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	cost(i,j) = ((branches_i[i].first != branches_j[j].first) +
		     (branches_i[i].second != branches_j[j].second));
    Alignement * align = EditDistance::applyMunkres(cost);
    double val =  diff + align->cost;
    delete align;
    return val;
  }else if(treelet_type == 9){
    int diff = code_i[0] != code_j[0];//v0
    diff +=  code_i[2] != code_j[2];//e0
    diff +=  code_i[4] != code_j[4];//v1
    //Permutations between e1v2e2v3 and e3v4e4v5
    //e1v2e2v3 -> e1v2e2v3, e3v4e4v5->e3v4e4v5
    int diff_1=0;
    for(int i=0;i<8;i++)
      diff_1 += (code_i[6+(i*2)] != code_j[6+(i*2)]);
    //e1v2e2v3 -> e3v4e4v5, e3v4e4v5->e1v2e2v3
    int diff_2 = (code_i[6] != code_j[14]);
    diff_2 += (code_i[8] != code_j[16]);
    diff_2 += (code_i[10] != code_j[18]);
    diff_2 += (code_i[12] != code_j[20]);
    diff_2 += (code_i[14] != code_j[6]);
    diff_2 += (code_i[16] != code_j[8]);
    diff_2 += (code_i[18] != code_j[10]);
    diff_2 += (code_i[20] != code_j[12]);

    return diff + (diff_1 < diff_2)?diff_1:diff_2;
    
  }else if(treelet_type == 11){
    int diff = 0;
    for(int i = 0; i<code_i.length()-(6*2);i++) //permutations between e2v3,e3v4 and e4v5
      diff += code_i[i] != code_j[i];
    vector<pair<char,char> > branches_i(3);
    vector<pair<char,char> > branches_j(3);
    branches_i[1] = pair<char,char>(code_i[code_i.length()-11],code_i[code_i.length()-9]);
    branches_i[0] = pair<char,char>(code_i[code_i.length()-7],code_i[code_i.length()-5]);
    branches_i[2] = pair<char,char>(code_i[code_i.length()-3],code_i[code_i.length()-1]);
    branches_j[1] = pair<char,char>(code_j[code_j.length()-11],code_j[code_j.length()-9]);
    branches_j[0] = pair<char,char>(code_j[code_j.length()-7],code_j[code_j.length()-5]);
    branches_j[2] = pair<char,char>(code_j[code_j.length()-3],code_j[code_j.length()-1]);
    CImg<double> cost(3,3);
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	cost(i,j) = ((branches_i[i].first != branches_j[j].first) +
		     (branches_i[i].second != branches_j[j].second));
    Alignement * align = EditDistance::applyMunkres(cost);
    return diff + align->cost;

  }else if(treelet_type == 12){
    int permut_intra_g12[8][21] = {
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}, //v0v1v2v3v4v5
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,18,15,20,17,14,19,16}, //v0v1v2v3v5v4
      {0,1,6,3,8,5,2,7,4,9,10,11,12,13,14,15,16,17,18,19,20}, //v0v2v1v3v4v5
      {0,1,6,3,8,5,2,7,4,9,10,11,12,13,18,15,20,17,14,19,16}, //v0v1v2v3v4v5
      {12,1,14,3,16,5,18,7,20,9,10,11,0,13,2,15,4,17,6,19,8}, //v3v4v5v0v1v2
      {12,1,14,3,16,5,18,7,20,9,10,11,0,13,6,15,8,17,2,19,4}, //v3v4v5v0v2v1
      {12,1,18,3,20,5,14,7,16,9,10,11,0,13,2,15,4,17,6,19,8}, //v3v5v4v0v1v2
      {12,1,18,3,20,5,14,7,16,9,10,11,0,13,6,15,8,17,2,19,4} //v3v5v4v0v2v1
    };
    int diff = INT_MAX;
    // string code_i_perms[8];
    for(int i=0;i<8;i++){
      // code_i_perms[i] = string(code_i);
      int cur_diff = 0;
      for(int j=0;j<21;j++)
	cur_diff += (code_i[permut_intra_g12[i][j]] != code_j[j]);
      diff = (diff<cur_diff)?diff:cur_diff;
    }
    return diff;
  }
  return 0;
  
}

double TreeletEditDistance::operator() (int t_i, std::string code_i, int t_j, std::string code_j){
  
  int k = mcs[t_i][t_j];
  vector<string> codes_i = Codes(t_i,code_i,k);
  vector<string> codes_j = Codes(t_j,code_j,k);
  double c_struct = edit_distance_struct[t_i][t_j];
  double c_subs = DBL_MAX;
  for(int i=0;i<codes_i.size();i++)
    for(int j=0;j<codes_j.size();j++)
      {
	double tmp = compareNodesCodes(codes_i[i],codes_j[j],k);
	c_subs = (c_subs<tmp)?c_subs:tmp;
      }
  return c_subs*substitution_cost + c_struct * deletion_cost;
}

double TreeletEditDistance::operator() ( pandore::Collection * c1, pandore::Collection * c2 ){
  return 0.0;
}


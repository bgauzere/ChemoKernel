/*
 * @file InterTreeletKernel.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Fri Apr  6 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "InterTreeletKernel.h"
#include "MoleculeGraph.h"
#include "TreeletEnumerator.h"
using namespace std;
using namespace cimg_library;
using namespace pandore;
#define EPS 0.0001

void InterTreeletKernel::regul(){
  CImg<double> eigvals;
  CImg<double> eigvects;
  this->sim_treelets.symmetric_eigen(eigvals,eigvects);
  cout << "regul ..." << endl;
  // for(int i=0;i<nb_treelets;i++)
  //   for(int j=0;j<nb_treelets;j++)
  //     assert(sim_treelets(i,j) == sim_treelets(j,i));
  if(eigvals[nb_treelets-1] < 0)
    for (unsigned int i=0; i<nb_treelets; ++i)
      sim_treelets(i,i) -= eigvals[nb_treelets-1];
  
}

/*XXX: Constructeur à factoriser*/
InterTreeletKernel::InterTreeletKernel(GraphKernel * k_inter_treelets, 
				       Kernel * k_inter_graphs, 
				       treelet_spectrum ** treelets, int nb_treelets){
  this->nb_treelets = nb_treelets;
  this->sim_treelets =  CImg<double>(nb_treelets,nb_treelets);
  this->k_inter_graphs = k_inter_graphs;
  Collection ** col_treelets = new Collection*[nb_treelets];
  int cur_treelet = 0;
  for(int i = 0; i < SIZE_SPECTRUM; i++) //Parcours des 13 Treelets
    {
      treelet_spectrum::iterator it = (*treelets)[i].begin();
      for(;it != (*treelets)[i].end();it++,cur_treelet++){
	//cout << i << "," << it->first << endl;
	col_treelets[cur_treelet] = TreeletEnumerator::TreeletToCollection(i,it->first);
      }
    }
  if(k_inter_treelets != NULL){
    for (int i=0; i < nb_treelets; i++)
      for (int j=i; j < nb_treelets; j++){
	sim_treelets(i,j) = sim_treelets(j,i) = (*k_inter_treelets)(col_treelets[i],col_treelets[j]); 
	if(sim_treelets(i,j) > EPS)
	  selected_pairs.push_back(pair<int,int>(i,j));
      }
    regul();
  }
  //sim_treelets.display();
}
  // cout << "Valeurs propres inter treeleets :" << endl;
  // CImg<double> eigvals;
  // CImg<double> eigvects;
  // sim_treelets.symmetric_eigen(eigvals,eigvects);
  // bool isPD = true;
  // for (unsigned int i=0; i<nb_treelets; ++i)
  //   {
  //     cout << eigvals[i] << endl;
  //     // cout << eigvects[i] << endl;
  //     //   if(eigvals[i] < 0)
  //     // 	return false;
  //     // 
  //     if(eigvals[i] < 0)
  // 	isPD = false;
  //   }
  // cout << "Fin Valeurs propres inter treeleets :" << endl;



InterTreeletKernel::InterTreeletKernel(GraphKernel * k_inter_treelets, 
				       Kernel * k_inter_graphs, 
				       Collection ** treelets, int nb_treelets){
  this->nb_treelets = nb_treelets;
  this->sim_treelets =  CImg<double>(nb_treelets,nb_treelets);
  cout << "Sim Treelets matrix allocated"<< endl;
  this->k_inter_graphs = k_inter_graphs;
  Collection ** col_treelets = new Collection*[nb_treelets];
  int cur_treelet = 0;
  for(int i = 0; i < nb_treelets; i++)
    {
      col_treelets[i] = treelets[i];
    }
  if(k_inter_treelets != NULL){
    for (int i=0; i < nb_treelets; i++)
      for (int j=i; j < nb_treelets; j++){
	sim_treelets(i,j) = sim_treelets(j,i) = (*k_inter_treelets)(col_treelets[i],col_treelets[j]);
	if(sim_treelets(i,j) > EPS)
	  selected_pairs.push_back(pair<int,int>(i,j));
      }
    regul();
  }
  // sim_treelets.display();
}

double InterTreeletKernel::operator()(Collection* c1, Collection* c2){
  double scalar = 0.0;
  //calcul matriciel
  CImg<double> * vec_1 = (CImg<double> * ) c1->GETVALUE("labeled_special",Long);
  CImg<double> * vec_2 = (CImg<double> * ) c2->GETVALUE("labeled_special",Long);
  //return (vec_1.transpose()*sim_treelets*vec_2).sum;
  double sum = 0.0;
  for(int n =0;n<selected_pairs.size();n++){
    int i = selected_pairs[n].first;
    int j = selected_pairs[n].second;
    sum += (*k_inter_graphs)((*vec_1)(i), (*vec_2)(j))*sim_treelets(j,i);
    if(i!= j)
      sum += (*k_inter_graphs)((*vec_1)(j), (*vec_2)(i))*sim_treelets(j,i);
  }
  // double sum2 =0.0;
  // for(int i=0;i<nb_treelets;i++)
  //   for(int j=0;j<nb_treelets;j++)
  //     sum2 += ((*k_inter_graphs)((*vec_2)(j), (*vec_1)(i)))*sim_treelets(i,j);
  // cout << "sum: "  << sum << endl;
  // cout << "sum2: "  << sum2 << endl;
  // cout << "diff: "  << sum - sum2 << endl;

  return sum;
}

void InterTreeletKernel::setSimTreelets(CImg<double> sim){
  this->sim_treelets = CImg<double>(sim);
  for (int i=0; i < nb_treelets; i++)
    for (int j=i; j < nb_treelets; j++){
      if(sim_treelets(i,j) > EPS)
	selected_pairs.push_back(pair<int,int>(i,j));
    }
}

CImg<double> * InterTreeletKernel::getSimTreelets(){
  return &sim_treelets;
}

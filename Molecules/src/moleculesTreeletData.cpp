/*
 * @file moleculesTreeletData.cpp
 *
 * This program computes the correlation between the Graphlet Distribution and the property of the Molecule.
 *
 * @author Benoit GAUZERE 
 *
 * @version 1.0.0 (2010-11-2)
 */
#include <climits>
#include <cfloat>
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include "MoleculeGraph.h"
#include <pandore.h>
#include "string_utils.h"
#include "GraphletCountKernel.h"
#include "TreeletEditDistance.h"
#include "MoleculeGraphEditDistanceMCS.h"
#include "MoleculesDataset.h"
#include "CImg.h"

using namespace std;
using namespace pandore;
using namespace cimg_library;

struct Options
{
  char* trainset_path;
  char* dataset_file;
  char* output_file;
  double cns;
  double cni;
  double sigma;
  bool display;
  int type_treelet;
  int nbAtomsMax;
  bool noUnique;
  bool noUseless;
  char * test_file;
  bool sim;
  bool approx;
  int labelCycles;
};

/**
 * Show the program's usage.
 */

void showUsage ()
{
  cerr << "Usage : moleculesTreeletData trainset_path dataset_file " << endl;
  cerr << " Options:" << endl;
  cerr << "\t -o file_prefix" << endl;
  cerr << "\t -sim : compute similarity between treelets" << endl;
  cerr << "\t -cni cni" << endl;
  cerr << "\t -cns cns" << endl;
  cerr << "\t -s sigma" << endl;
  cerr << "\t -a : Compute an approximate edit distance" << endl;

  cerr << "\t -d no output, display sub gram matrices" << endl;
  cerr << "\t -t Type_Treelet :" << endl;
  cerr << "\t \t 0 Treelet" << endl;
  cerr << "\t \t 1 Treelet d'hypergraphe de cycles contractés" << endl;
  cerr << "\t1 -- V1 for acyclic - Parcours par couche pour cas prenant en compte les cycles" << endl;
  cerr << "\t2 -- V2 for acyclic - Parcours par branche pour cas prenant en compte les cycles" << endl;
  cerr << "\t3 -- Parcours par couche pour cas prenant en compte les cycles + Extension sans treelets" << endl;
  cerr << "\t4 -- Parcours par branche pour cas prenant en compte les cycles + Extension sans treelets" << endl;
  cerr << "\t -nb nombre maximum d'atome par treelets" << endl;
  cerr << "\t -u supprimer les treelets présents dans une seul molécule" << endl;
  cerr << "\t -p supprimer les treelets présents dans toutes les molécule le même nombre de fois" << endl;
  cerr << "\t -test dataset_test" << endl;
  cerr << "-lC Modifications of label"<< endl;
  cerr << "\t0 -- No modifications" << endl;
  cerr << "\t1 -- Number of relevant cycles add to the label of an atom" << endl;
  cerr << "\t2 -- Size of each relevant cycles add to the label of an atom" << endl;
}

void readOptions (int argc, char** argv, Options* options)
{
  for(int i = 0; i < argc; i++){
    cout << argv[i] << " ";
  }
  cout << endl;
  if (argc < 3)
    {
      showUsage();
      exit(1);
    }
  options->trainset_path = argv[1];
  options->dataset_file = argv[2];
  options->display = false;
  options->type_treelet = 0;
  options->nbAtomsMax=SIZE_MAX;
  options->noUnique=false;
  options->noUseless=false;
  options->test_file=NULL;
  options->sim=false;
  options->approx= false;
  options->labelCycles=0;
  int i=3;
  while (i<argc)
    {
      if (strncmp(argv[i], "-o", 2) == 0)
	{
	  options->output_file = argv[i+1];
	  i+=2;
	}
      else if (strcmp(argv[i], "-test") == 0)
	{
	  options->test_file = argv[i+1];
	  i+=2;
	}
      else if (strcmp(argv[i], "-sim") == 0)
	{
	  options->sim=true;
	  i+=1;
	}
      else if (strcmp(argv[i], "-a") == 0)
	{
	  options->approx=true;
	  i+=1;
	}
      else if (strcmp(argv[i], "-cni") == 0)
	{
	  options->cni = atof(argv[i+1]);
	  i+=2;
	}
      else if (strcmp(argv[i], "-cns") == 0)
	{
	  options->cns = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-s", 2) == 0)
	{
	  options->sigma = atof(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-t", 2) == 0)
	{
	  options->type_treelet = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-d", 2) == 0)
	{
	  options->display = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-nb", 2) == 0)
	{
	  options->nbAtomsMax = atoi(argv[i+1]);
	  i+=2;
	}
      else if (strncmp(argv[i], "-u", 2) == 0)
	{
	  options->noUnique = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-p", 2) == 0)
	{
	  options->noUseless = true;
	  i+=1;
	}
      else if (strncmp(argv[i], "-lC", 3) == 0)
	{
	  options->labelCycles = atoi(argv[i+1]);
	  i+=2;
	} 
      else
	{
	  cerr << "Unknown parameter : " << argv[i] << endl;
	  showUsage();
	  exit(1);
	}
    }
}


bool comp (string l, string r)
{ 
  return (l.compare(r) < 0); 
}

bool param_condition(double x)
{
  return true;
}

int main (int argc, char** argv)
{	
  Options options;
  readOptions(argc, argv, &options);

  MoleculeGraph::initTable();
  MoleculesDataset dataset (options.trainset_path, options.dataset_file);
  MoleculesDataset testset;
  
  if(options.test_file)
    testset =  MoleculesDataset(options.trainset_path, options.test_file);
	
  /**** Labeled Treelet Distribution ****/
  if(options.type_treelet ==0){
    long diff;	
    clock_t start,end;
    struct timeval s_start, s_end;
    start = gettimeofday(&s_start,NULL);
    dataset.computeLabeledSpectrums();
    end = gettimeofday(&s_end,NULL);//times(&s_end);
    diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
      (1000000*s_start.tv_sec + s_start.tv_usec);
    cout.precision(DBL_DIG);
    cout << "Treelet computed :" << diff << endl;

    if(options.test_file)
      testset.computeLabeledSpectrums();
  }
  
  else if (options.type_treelet ==1){
    dataset.computeCCHypergraph();
    dataset.computeContractedCycleSpectrums();
    if(options.test_file){
      testset.computeCCHypergraph();
      testset.computeContractedCycleSpectrums();
    }
  }
  else{
    cerr << "Type treelet non défini"<< endl;
    return -1;
  }
    
  
  /**** Init Files Operations ****/
  char * treelets_codes_suffix = "_treelets_codes";
  char * treelets_codes_filename = new char[strlen(options.output_file)+strlen(treelets_codes_suffix)+1];
  treelets_codes_filename = strcpy(treelets_codes_filename,options.output_file);
  treelets_codes_filename = strcat(treelets_codes_filename,treelets_codes_suffix);
  ofstream treelets_codes_outfile (treelets_codes_filename);	

  char * treelets_distrib_suffix = "_treelets_distrib";
  char * treelets_distrib_filename = new char[strlen(options.output_file)+strlen(treelets_distrib_suffix)+1];
  treelets_distrib_filename = strcpy(treelets_distrib_filename,options.output_file);
  treelets_distrib_filename = strcat(treelets_distrib_filename,treelets_distrib_suffix);
  ofstream treelets_distrib_outfile (treelets_distrib_filename);

  /* Fichier qui garde les types des treelets
   * Il est de la forme: 
   * num du dernier treelet de type 0
   * num du dernier treelet de type 1
   * num du dernier treelet de type 2
   * ....
   */
  char * treelets_types_suffix = "_treelets_types";
  char * treelets_types_filename = new char[strlen(options.output_file)+strlen(treelets_types_suffix)+1];
  treelets_types_filename = strcpy(treelets_types_filename,options.output_file);
  treelets_types_filename = strcat(treelets_types_filename,treelets_types_suffix);
  ofstream treelets_types_outfile (treelets_types_filename);	
	

  ofstream * treelets_test_distrib_outfile;
  ofstream * y_test_outfile;
  if(options.test_file){
    char * treelets_test_distrib_suffix = "_treelets_test_distrib";
    char * treelets_test_distrib_filename = new char[strlen(options.output_file)+strlen(treelets_test_distrib_suffix)+1];
    treelets_test_distrib_filename = strcpy(treelets_test_distrib_filename,options.output_file);
    treelets_test_distrib_filename = strcat(treelets_test_distrib_filename,treelets_test_distrib_suffix);
    treelets_test_distrib_outfile =  new ofstream(treelets_test_distrib_filename);	
    
    char * y_test_suffix = "_y_test";
    char * y_test_filename = new char[strlen(options.output_file)+strlen(y_test_suffix)+1];
    y_test_filename = strcpy(y_test_filename,options.output_file);
    y_test_filename = strcat(y_test_filename,y_test_suffix);
    y_test_outfile = new ofstream(y_test_filename);	
  }
  

  char * y_suffix = "_y";
  char * y_filename = new char[strlen(options.output_file)+strlen(y_suffix)+1];
  y_filename = strcpy(y_filename,options.output_file);
  y_filename = strcat(y_filename,y_suffix);
  ofstream y_outfile (y_filename);

  //Initialisation de l ensemble des treelets présents
  int nb_treelets = 0;
  int nb_type = SIZE_SPECTRUM;
  if(options.type_treelet ==4)
    nb_type=SIZE_SPECTRUM*2; //On discerne les codes des treelets et des stereo, afin de pouvoir plus efficacement travailler dessus en matlab
  map<string,int, bool(*)(string,string)> * treelets = new map<string,int, bool(*)(string,string)>[nb_type];
  for(int j=0;j<nb_type;++j)
    treelets[j] = map<string,int, bool(*)(string,string)>(comp);
  treelet_spectrum * distribution;
  treelet_spectrum * distribution2;
  if(options.type_treelet ==0)
    distribution = dataset.getTreeletDistribution(param_condition);
  else if(options.type_treelet ==1)
    distribution = dataset.getContractedCycleTreeletDistribution(param_condition);
  
  for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
    {
      treelet_spectrum::iterator it = distribution[j].begin();
      for(;it != distribution[j].end();it ++)
	{
	  // char O_char[2];
	  // O_char[0]=8;
	  // O_char[1]='\0';
	  // if( strstr((it->first).c_str(),O_char)!=NULL)
	    if(treelets[j].find(it->first) == treelets[j].end())
	      {
		treelets[j].insert(pair<string,int>(it->first,0)); 
		nb_treelets ++;
	      }
	}
    }

  if(options.noUseless)
    for(int j=0;j<SIZE_SPECTRUM;++j)
      {
	map<string,int, bool(*)(string,string)>::iterator it = treelets[j].begin();
	while (it != treelets[j].end())
	  {
	    bool to_suppr=true;
	    int nbOccur=0;
	    treelet_spectrum ** spectrum;
	    if(options.type_treelet == 0)
	      spectrum = dataset.getTreeletSpectrum(0);
	    else if(options.type_treelet == 1)
	      spectrum = dataset.getContractedCycleSpectrum(0);
	    
	    if(spectrum[j]->find(it->first)!=spectrum[j]->end())
	      {
		nbOccur=(spectrum[j]->find(it->first))->second;
	      }
	    else
	      {
		to_suppr=false;
	      }

	    for(unsigned int i=1;i<dataset.size() && to_suppr;i++)
	      {
		if(options.type_treelet == 0)
		  spectrum = dataset.getTreeletSpectrum(i);
		else if(options.type_treelet == 1)
		  spectrum = dataset.getContractedCycleSpectrum(i);
		
		if(spectrum[j]->find(it->first)!=spectrum[j]->end())
		  {
		    if(nbOccur!=(spectrum[j]->find(it->first))->second)
		      to_suppr=false;
		  }
		else
		  {
		    to_suppr=false;
		  }
	      }
	    it++;
	    if(to_suppr)
	      {
		map<string,int, bool(*)(string,string)>::iterator copie_it=it;
		copie_it--;
		treelets[j].erase(copie_it);
	      }
	  }
      }

  char HO[6];


  
  //On associe chaque code a un index
  int * type_treelet = new int[nb_treelets];// Pour TreeletEditDistance
  string* code_treelet = new string[nb_treelets];// Pour TreeletEditDistance
  int index = 0;
  for(int j=0;j<nb_type;++j)
    {
      map<string,int, bool(*)(string,string)>::iterator it = treelets[j].begin();
      for(;it != treelets[j].end();it ++)
	{
	  stringstream key;
	  key <<  j << "-" << it->first;
	  treelets_codes_outfile << key.str() << endl;
	  it->second = index;
	  type_treelet[index] = j;
	  code_treelet[index] = it->first;
	  index ++;
	}
      treelets_types_outfile << index << endl;
    }
  
  vector< vector<int> > molecules;// =  vector< vector<int> >(dataset.size());
  for(unsigned int i=0;i<dataset.size();i++)
    {
      vector<int> current_molecule = vector<int>(nb_treelets,0);
      treelet_spectrum ** spectrum;
      treelet_spectrum ** spectrum2;
      if(options.type_treelet == 0)
	spectrum = dataset.getTreeletSpectrum(i);
      else if(options.type_treelet == 1)
	spectrum = dataset.getContractedCycleSpectrum(i);
      	
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	{
	  treelet_spectrum::iterator it = spectrum[j]->begin();
	  for(;it != spectrum[j]->end();it ++)
	    {
	      if(treelets[j].find(it->first)!=treelets[j].end()){
		int index = treelets[j].find(it->first)->second;
		current_molecule[index] = it->second;
	      }
	    }
	}
      if(options.type_treelet == 4)
	for(int j=0;j<SIZE_SPECTRUM;++j)
	  {
	    treelet_spectrum::iterator it = spectrum2[j]->begin();
	    for(;it != spectrum2[j]->end();it ++)
	      {
		if(treelets[j+SIZE_SPECTRUM].find(it->first)!=treelets[j+SIZE_SPECTRUM].end()){
		  int index = treelets[j+SIZE_SPECTRUM].find(it->first)->second;
		  current_molecule[index] = it->second;
		}
	      }
	  }

      molecules.push_back(vector<int>(current_molecule));
      for(unsigned int j=0;j<current_molecule.size();j++){
	treelets_distrib_outfile  << current_molecule[j] ;
	if (j != current_molecule.size() -1)
	  treelets_distrib_outfile<< ",";
	else
	  treelets_distrib_outfile<< endl;
      }
      y_outfile << dataset.getParameter(i) << endl;
    }
  treelets_distrib_outfile.close();
  treelets_codes_outfile.close();
  treelets_types_outfile.close();
  y_outfile.close();
  cout << "Distribution treelets apprentissage ok  ..." << nb_treelets << "treelets"<< endl;
  if(options.test_file != NULL){
    vector< vector<int> > molecules_test;// =  vector< vector<int> >(dataset.size());
    for(unsigned int i=0;i<testset.size();i++)
      {
	vector<int> current_molecule = vector<int>(nb_treelets,0); // on remplit le spectrum par rapport au trainset
	treelet_spectrum ** spectrum;
	if(options.type_treelet == 0)
	  spectrum = testset.getTreeletSpectrum(i);
	else if(options.type_treelet == 1)
	  spectrum = testset.getContractedCycleSpectrum(i);
	
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  {
	    treelet_spectrum::iterator it = spectrum[j]->begin();
	    for(;it != spectrum[j]->end();it ++)
	      {
		if(treelets[j].find(it->first)!=treelets[j].end()){
		  int index = treelets[j].find(it->first)->second;
		  current_molecule[index] = it->second;
		}
	      }
	  }
	molecules_test.push_back(vector<int>(current_molecule));
	// // cout << molecules_test.size() << endl;
	for(unsigned int j=0;j<current_molecule.size();j++){
	  (*treelets_test_distrib_outfile)  << current_molecule[j] ;
	  if (j != current_molecule.size() -1)
	    (*treelets_test_distrib_outfile)<< ",";
	  else
	    (*treelets_test_distrib_outfile)<< endl;
	}
	(*y_test_outfile) << testset.getParameter(i) << endl;
      }
    treelets_test_distrib_outfile->close();
    y_test_outfile->close();
    cout << "Distribution treelets test ok  ..." << nb_treelets << "treelets"<< endl;
  }



  if(options.sim){
    /*Similarité entre treelets que pour treelets classiques*/
    if(options.type_treelet == 0) { 
      char * treelets_sim_suffix = "_treelets_sim";
      char * treelets_sim_filename = new char[strlen(options.output_file)+strlen(treelets_sim_suffix)+1];
      treelets_sim_filename = strcpy(treelets_sim_filename,options.output_file);
      treelets_sim_filename = strcat(treelets_sim_filename,treelets_sim_suffix);
      ofstream treelets_sim_outfile (treelets_sim_filename);	
      treelets_sim_outfile.precision(20);

      char * treelets_distance_suffix = "_treelets_distance";
      char * treelets_distance_filename = new char[strlen(options.output_file)+strlen(treelets_distance_suffix)+1];
      treelets_distance_filename = strcpy(treelets_distance_filename,options.output_file);
      treelets_distance_filename = strcat(treelets_distance_filename,treelets_distance_suffix);
      ofstream treelets_distance_outfile (treelets_distance_filename);	
      treelets_distance_outfile.precision(20);


      Collection ** col_treelets = new Collection*[nb_treelets];      
      for(int i = 0; i < nb_treelets; i++) //Parcours des 13 Treelets
	col_treelets[i] = TreeletEnumerator::TreeletToCollection(type_treelet[i],code_treelet[i]);
      
      TreeletEditDistance::Init();

      clock_t itk_start,itk_end;
      struct timeval itk_s_start, itk_s_end;
      itk_start = gettimeofday(&itk_s_start,NULL);
      EditDistance * ed;
      if (options.approx)
	ed = new MoleculeGraphEditDistanceMCS(options.cns, options.cni);
      else
	ed = new TreeletEditDistance(options.cns,options.cni);
      
      TreeletEditDistance ted(options.cns,options.cni);
      
      TreeletEditDistance::Init();
      
      CImg<double> D_treelets(nb_treelets,nb_treelets,1,1,0.0);
      CImg<double> sim_treelets(nb_treelets,nb_treelets,1,1,0.0);
      for(int i=0;i<nb_treelets;i++)
	for(int j=i+1;j<nb_treelets;j++){
	  if(!options.approx)
	    D_treelets(i,j) = D_treelets(j,i)  = ted(type_treelet[i], code_treelet[i],
						     type_treelet[j], code_treelet[j]);
	  else
	    D_treelets(i,j) = D_treelets(j,i)  = (*ed)(col_treelets[i],col_treelets[j]);
	  
	  sim_treelets(j,i) =sim_treelets(i,j)=  exp(- pow(D_treelets(i,j),2)/(options.sigma));
	  // if(sim_treelets(i,j) < 0.001)
	  if(sim_treelets(i,j) < DBL_EPSILON)
	    sim_treelets(i,j) = 0.0;
	  //cout << sim_treelets(i,j);
	  // if(j!=nb_treelets-1)
	  //   cout <<  ",";
	  // else
	  //   cout << endl;
	}
      //	sim_treelets.display();
      //	Regul
      itk_end = gettimeofday(&itk_s_end,NULL);//times(&s_end);
      long itk_time = (1000000*itk_s_end.tv_sec + itk_s_end.tv_usec) - 
	(1000000*itk_s_start.tv_sec + itk_s_start.tv_usec);	

      cout << "Similarité inter treelets ok ..." << endl;
      cout << "Similarité inter treelets calculée en  " << itk_time << "micro seconds."<< endl;    

      // cout << "Regularisation ..." << endl;
      CImg<double> eigvals;
      CImg<double> eigvects;
      sim_treelets.symmetric_eigen(eigvals,eigvects);
      // for(int i=0;i<nb_treelets;i++)
      //   for(int j=0;j<nb_treelets;j++)
      //     assert(sim_treelets(i,j) == sim_treelets(j,i));
    
      if(eigvals[nb_treelets-1] < 0)
	for (unsigned int i=0; i<nb_treelets; ++i)
	  sim_treelets(i,i) -= eigvals[nb_treelets-1];
      // sim_treelets.symmetric_eigen(eigvals,eigvects);
      // for (unsigned int i=0; i<nb_treelets; ++i)
      // 	cout << eigvals[i] << endl;
      for(int i=0;i<nb_treelets;i++){
	for(int j=0;j<nb_treelets;j++){
	  treelets_sim_outfile << sim_treelets(i,j);
	  treelets_distance_outfile << D_treelets(i,j);		  
	  if(j!=nb_treelets-1){
	    treelets_distance_outfile << ",";
	    treelets_sim_outfile << ",";
	  }
	  else
	    {
	      treelets_distance_outfile << endl;
	      treelets_sim_outfile << endl;
	    }
	}
      }
    
      if(options.display){
	for(int i=0;i<nb_treelets;i++)
	  for(int j=0;j<nb_treelets;j++){
	    //paire i,j
	    cout << i << "," << j << endl;
	    sim_treelets = CImg<double> (molecules.size(),molecules.size());
	    for(int n=0;n<molecules.size();n++){
	      for(int m=0;m<molecules.size();m++){
		sim_treelets(n,m) = 
		  exp(- pow(molecules[n][i]-molecules[m][j],2)/(3*3)) +
		  exp(- pow(molecules[n][j]-molecules[m][i],2)/(3*3));
		exp(- pow(molecules[n][i]-molecules[m][i],2)/(3*3)) +
		  exp(- pow(molecules[n][j]-molecules[m][j],2)/(3*3));
		if(n==m)
		  cout << sim_treelets(n,m) << "," << endl;
	      }
	    }
	  
	    sim_treelets.display();
	  }
      }
    }
  }
  dataset.destroy();
  if(options.test_file)
    testset.destroy();
  
  return 0;
}

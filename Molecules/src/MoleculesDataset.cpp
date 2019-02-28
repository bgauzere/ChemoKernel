/*
 * @file MoleculeGraph.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */


#include <iostream>
#include <fstream>
#include <time.h>
#include <assert.h>
#include <float.h>
#include <map>
#include <libgen.h>

#include "MoleculesDataset.h"
#include "string_utils.h"

using namespace std;


MoleculesDataset::MoleculesDataset(const MoleculesDataset& d) : Dataset(d)
{
  _molecules = d._molecules;
  _filenames = d._filenames;
  nbTreelets = d.nbTreelets;
}

MoleculesDataset::MoleculesDataset(const char* path, const char* filename)
{
  loadDataset(path, filename);
}

MoleculesDataset::MoleculesDataset(const char* filename)
{

  loadDataset(dirname((char*)filename), basename((char*)filename));
}

void MoleculesDataset::loadDataset (const char* path, const char* filename)
{
  char* dataset;
	
  if (path[strlen(path)-1] != '/')
    {
      dataset = new char[strlen(path)+strlen(filename)+2];
      sprintf(dataset, "%s/%s", path, filename);
    }
  else
    {
      dataset = new char[strlen(path)+strlen(filename)+1];
      sprintf(dataset, "%s%s", path, filename);
    }
	
  ifstream f (dataset);
  long id = 0;
  if (f)
    {
      string s;
		
      while (getline(f, s))
	{
	  if (s[0] != '#')
	    {
	      const char * s_c_str = s.c_str();
	      const char * pch=strrchr(s_c_str,' ');
	      char* ctfile = new char [pch-s_c_str+1];
	      ctfile = strncpy(ctfile,s_c_str,pch-s_c_str);
	      ctfile[pch-s_c_str] = '\0';
	      
	      char * parameter = new char[strlen(pch)+1];
	      strcpy(parameter,pch+1);
	      parameter[strlen(pch)] = '\0';
	      char* ctfile_absolute = new char [strlen(path)+strlen(ctfile)+2];
	      ctfile_absolute = strcpy(ctfile_absolute,path);
	      ctfile_absolute = strcat(ctfile_absolute,"/");
	      ctfile_absolute = strcat(ctfile_absolute,ctfile);
	      
	      // Creation of the MoleculeGraph
	      MoleculeGraph* tmp = new MoleculeGraph (ctfile_absolute,id);
	      delete [] ctfile_absolute;
	      if (tmp->getCollection() != NULL)//Molecule graph success
		{
		  id++;
		  _collections.push_back(tmp->getCollection());
		  _parameters.push_back(atof(parameter));
		  _molecules.push_back(tmp);
		  _filenames.push_back(ctfile); 	      
		}
	      delete [] parameter;
	    }
	}      
    }
  else
    cerr << "Error : can't load the dataset" << endl;
	
  f.close();
	
  delete[] dataset;
}

void MoleculesDataset::destroy ()
{
  for (unsigned int i=0; i<_molecules.size(); ++i)
    delete _molecules[i];
}

void MoleculesDataset::reduceToN(unsigned int new_size){
  while (_collections.size() > new_size){
    //
    srand(time(NULL));
    int k = rand() % _collections.size();
    _molecules.erase(_molecules.begin() + k);
    _collections.erase(_collections.begin() + k
);
    _parameters.erase(_parameters.begin() + k);
  }
  assert(this->size() == new_size);
}

void MoleculesDataset::eraseSome(list<int> to_erase){
  //sort List
  to_erase.sort();
  //erase higher to lower
  list<int>::iterator it = to_erase.end();
  do{
    it --;
    int k = *it;
    _molecules.erase(_molecules.begin() + k);
    _collections.erase(_collections.begin() + k);
    _parameters.erase(_parameters.begin() + k);
    _filenames.erase(_filenames.begin() + k);
  }while(it != to_erase.begin());
}

void MoleculesDataset:: eraseByFilename(std::vector<char *> to_erase)
{
  list<int> to_eraseList;
  for(vector<char *>::iterator it = to_erase.begin();it != to_erase.end(); ++it)
    {
      to_eraseList.push_back(this->getNumberByFilename(*it));
    }
  if(to_eraseList.size()!=0)
    this->eraseSome(to_eraseList);
}

void MoleculesDataset::computeSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeGraphletSpectrum();
  }  
}

void MoleculesDataset::computeLabeledSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeLabeledTreeletSpectrum();
  }  
}
void MoleculesDataset::computeCyclesSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeCycleTreeletSpectrum();
  }  
}

void MoleculesDataset::computeCCHypergraph(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeContractedCycleGraphs();
  }  
}
void MoleculesDataset::computeAugmentedCycles(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeAugmentedCycles();
  }  
}


void MoleculesDataset::computeContractedCycleSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeContractedCycleSpectrum();
  }  
}

void MoleculesDataset::computeAugmentedCycleSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    (*it)->computeAugmentedCycleSpectrum();
  }  
}

void MoleculesDataset::printSpectrums(){
  int n = 0;
  cout << "Graphlet :\t" << "G0\t" <<"G1\t" <<"G2\t" <<"G3\t" 
       <<"G4\t" <<"G5\t" <<"G6\t" <<"G7\t" <<"G8\t" <<"G9\t" <<"G10\t" << "G11\t" <<"G12\t" <<"G13\t" << endl;
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it){
    double * spectrum = (*it)->getSpectrum();
    cout << "Molecule " <<  n << " :\t";
    for(int i = 0; i < (*it)->getSpectrumSize();i++)
      cout << spectrum[i] << "\t" ;
    cout << endl;
    n++;
  }
}

void MoleculesDataset::printLabeledSpectrums(){
  int n =0;
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it)
    {
      treelet_spectrum ** spectrum = (*it)->getLabeledTreeletSpectrum();
      cout << "Molecule " <<  n << " :" << endl;
      //XXX: Virer ce 13
      for(int i = 0; i < 13;i++)
	{
	  cout << "Graphlet " << i << " : " << endl;
	  treelet_spectrum::iterator spec_it;
	  for (spec_it=spectrum[i]->begin() ; spec_it != spectrum[i]->end(); spec_it++)
	    cout << "\t" << (*spec_it).first << " : " << (*spec_it).second << endl;
	}
      n++;
    }
}

void MoleculesDataset::normalizeLabeledSpectrums(){
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();it != _molecules.end(); ++it)
    {
      (*it)->normalizeLabeledSpectrum();
    }
}


void MoleculesDataset::computeGraphletCorrelation(){
 
  _correlation = new double[SIZE_SPECTRUM];
  this->computeSpectrums();
  
  int N = this->size();
  
  double temp_moy = 0.0;
  double * freq_moy = new double[SIZE_SPECTRUM];
  memset(freq_moy,0,sizeof(double)*SIZE_SPECTRUM);
  
  //Calcul des moyennes
  for (int i = 0; i < N ; ++i)
    { 
      //Parcours des molecules du Dataset
      double * spectrum = _molecules[i]->getSpectrum();
      temp_moy += _parameters[i];
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque Graphlet
	freq_moy[j] += spectrum[j];
    }
  temp_moy /= (double)N;
  for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque Graphlet
    freq_moy[j] /= (double)N;
	
  double * correlation_num = new double[SIZE_SPECTRUM];
  memset(correlation_num,0,sizeof(double)*SIZE_SPECTRUM);

  double * correlation_den_x = new double[SIZE_SPECTRUM];
  memset(correlation_den_x,0,sizeof(double)*SIZE_SPECTRUM);

  double * correlation_den_y = new double[SIZE_SPECTRUM];
  memset(correlation_den_y,0,sizeof(double)*SIZE_SPECTRUM);

  for (int i=0;i<N;++i)
    { 
      double * spectrum = _molecules[i]->getSpectrum();
      double temp = _parameters[i];

      for(int j=0;j<SIZE_SPECTRUM;++j){//Parcours de chaque Treelet
	correlation_num[j] += (spectrum[j]-freq_moy[j])*(temp-temp_moy);
	correlation_den_x[j] += pow((spectrum[j]-freq_moy[j]),2);
	correlation_den_y[j] += pow((temp-temp_moy),2);
      }
    }
	
  for(int j=0;j<SIZE_SPECTRUM;++j)
    {//Parcours de chaque Graphlet
      _correlation[j] = correlation_num[j]/(sqrt(correlation_den_x[j])*sqrt(correlation_den_y[j]));
    }
  
}

void MoleculesDataset::computeLabeledGraphletCorrelation(){
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;

  this->computeLabeledSpectrums();

  double N = this->size();
  double temp_moy = 0.0;
  vector<double> * mean_freq = new vector<double>[SIZE_SPECTRUM];
  //Associe un code de treelet a un indice dans mean_freq
  treelet_spectrum * code_to_indice = 
    new treelet_spectrum[SIZE_SPECTRUM]; //v^t
  for(int i=0;i<SIZE_SPECTRUM;i++)
    {
      mean_freq[i] = vector<double>();
      code_to_indice[i] = treelet_spectrum(fn_pt);
   }
  
  //Calcul des moyennes
  int nb_treelet_type = 0;
  for (int i = 0; i < N ; ++i) //Parcours des molécules
    { 
      temp_moy += _parameters[i];
      
      treelet_spectrum ** spectrum = _molecules[i]->getLabeledTreeletSpectrum();
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de Treelet
	for(treelet_spectrum::iterator it = spectrum[j]->begin();
	    it!=spectrum[j]->end(); it ++) //Parcours de chaque Treelet
	  {
	    int indice;
	    //Creation ou recuperation de l'indice dans code_to_indice
	    if(code_to_indice[j].find(it->first) != code_to_indice[j].end())
	      indice = code_to_indice[j].find(it->first)->second;
	    else
	      {
		mean_freq[j].push_back(0.0);
		nb_treelet_type++;
		code_to_indice[j].insert(pair<string, int>(it->first,indice=mean_freq[j].size()-1));
	      }
	    mean_freq[j][indice] += it->second;
	  }
    }

  temp_moy /= (double)N;
  for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
    for(unsigned int k=0; k<mean_freq[j].size();k++)//Parcours de chaque treelet
      mean_freq[j][k] /= (double)N;

  vector<double> * correlation_num = new vector<double>[SIZE_SPECTRUM];
  vector<double> * correlation_den_x = new vector<double>[SIZE_SPECTRUM];
  vector<double> * correlation_den_y = new vector<double>[SIZE_SPECTRUM];
  for(int i=0;i<SIZE_SPECTRUM;i++)
    {
      correlation_num[i] = vector<double>(mean_freq[i].size());
      correlation_den_x[i] = vector<double>(mean_freq[i].size());
      correlation_den_y[i] = vector<double>(mean_freq[i].size());
    }
  
  for (int i=0;i<N;++i)//Parcours de chaque molecule
    { 
      treelet_spectrum ** spectrum = _molecules[i]->getLabeledTreeletSpectrum();
      double temp = _parameters[i];
      for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de Treelet
	{
	  treelet_spectrum::iterator it_treelet;
	  for(it_treelet = code_to_indice[j].begin();
	      it_treelet != code_to_indice[j].end();
	      it_treelet++)
	    {//Parcours de chaque Treelet
	      double x_i = 0.0;
	      if(spectrum[j]->find(it_treelet->first) != spectrum[j]->end())
		x_i = spectrum[j]->find(it_treelet->first)->second;
	      int k = it_treelet->second;
	      correlation_num[j][k] += (x_i - mean_freq[j][k])*(temp-temp_moy);
	      correlation_den_x[j][k] += pow((x_i - mean_freq[j][k]),2);
	      //Inutile, commun a tout le monde
	      correlation_den_y[j][k] += pow(temp-temp_moy,2);
	    }
	}
    }

  _labeledCorrelation = new treelet_spectrum[SIZE_SPECTRUM];
  for(int j=0;j<SIZE_SPECTRUM;++j)
    {//Parcours de chaque type de Treelet
      _labeledCorrelation[j] = treelet_spectrum(fn_pt);
      treelet_spectrum::iterator it_treelet;
      for(it_treelet = code_to_indice[j].begin();
	  it_treelet != code_to_indice[j].end();
	  it_treelet++)
	{
	  double num = correlation_num[j][it_treelet->second]/N;
	  double den_x = correlation_den_x[j][it_treelet->second];
	  double den_y = correlation_den_y[j][it_treelet->second];
	  _labeledCorrelation[j].insert(pair<string, double>(it_treelet->first,
							     num/(sqrt(den_x/N)*sqrt(den_y/N))));
	}
    }  
}

treelet_spectrum * MoleculesDataset::getLabeledCorrelation(){
  return _labeledCorrelation;
}

double MoleculesDataset::getCorrelationCoeff(int graphlet){
  return _correlation[graphlet];
}

vector<MoleculesDataset*> MoleculesDataset::splitBySize(){
  //_collections, _parameters et _molecules a remplir pour chaque Dataset;
  
  //Recherche Taille max
  int size_max = 0;
  for(deque<MoleculeGraph*>::iterator it =  _molecules.begin();it != _molecules.end();it++)
    size_max = ((*it)->nbAtoms()>size_max)?(*it)->nbAtoms():size_max;
  
    

  
  
  vector<MoleculesDataset*> datasets(size_max+1); //datasets[0] non utilisé
  for(int i = 0; i <= size_max;i++)
    datasets[i] = new MoleculesDataset();

  int i = 0;
  for(deque<MoleculeGraph *>::iterator it = _molecules.begin();
      it != _molecules.end(); it ++)
    {
      int classe =  (*it)->nbAtoms();
      if(classe < 7)
	classe = 0;
      datasets[classe]->simpleAdd(_collections[i],
				  _molecules[i],
				  _parameters[i]);
      i++;
    }
  return datasets;
}


void MoleculesDataset::simpleAdd(pandore::Collection* c, MoleculeGraph * m ,double parameter){
  this->_collections.push_back(c);
  this->_molecules.push_back(m);
  this->_parameters.push_back(parameter);
}

void MoleculesDataset::normalizeParams(double fact){
  if(fact != 0)
    {
      for (unsigned int i = 0; i<_parameters.size();i++)
	_parameters[i] /= fact;    
    }
  else
    {
      //Normalisation entre 0 et 1
      double min_p = DBL_MAX;
      double max_p = DBL_MIN;
      
      for (unsigned int i = 0; i<_parameters.size();i++)
	{
	  min_p = min(_parameters[i],min_p);
	  max_p = max(_parameters[i],max_p);
	}
      for (unsigned int i = 0; i<_parameters.size();i++)
	{
	  _parameters[i] += (-min_p);
	  _parameters[i] /= (-min_p + max_p);
	}
    }
}


double * MoleculesDataset::getUnlabeledTreeletDistribution(bool (*param_condition)(double i)){
  double * distribution = new double[SIZE_SPECTRUM];
  memset(distribution,0,SIZE_SPECTRUM*sizeof(double));

  int N = _molecules.size();
  for(int i=0;i<N;i++)
    {
      if(param_condition(_parameters[i]))
	{
	  double * mol_distrib =_molecules[i]->getSpectrum();
	  for(int j=0;j<SIZE_SPECTRUM;j++)
	    {
	      distribution[j] += mol_distrib[j];
	    }
	}
    }
  return distribution;
}


void MoleculesDataset::computeRelevantCycles(){
   int N = _molecules.size();
   for(int i=0;i<N;i++)
     _molecules[i]->computeRelevantCycles();
}

void MoleculesDataset::computeSimpleCycles(int k){
   int N = _molecules.size();
   for(int i=0;i<N;i++)
     _molecules[i]->computeSimpleCycles(k);
}

vector<string> MoleculesDataset::getSimpleCycleCodes(int molecule){
  return _molecules[molecule]->getSimpleCycleCodes();
}


treelet_spectrum ** MoleculesDataset::getTreeletSpectrum(int molecule){
  return _molecules[molecule]->getLabeledTreeletSpectrum();
}
treelet_spectrum ** MoleculesDataset::getCycleTreeletSpectrum(int molecule){
  return _molecules[molecule]->getCycleTreeletSpectrum();
}

treelet_spectrum ** MoleculesDataset::getContractedCycleSpectrum(int molecule){
  return _molecules[molecule]->getContractedCycleSpectrum();
}
treelet_spectrum ** MoleculesDataset::getAugmentedCycleSpectrum(int molecule){
  return _molecules[molecule]->getAugmentedCycleSpectrum();
}

treelet_spectrum * MoleculesDataset::getTreeletDistribution(bool (*param_condition)(double i)){
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  nbTreelets = 0;
  treelet_spectrum * distribution = new treelet_spectrum[SIZE_SPECTRUM];
  for(int i=0;i<SIZE_SPECTRUM;i++)
    distribution[i] = treelet_spectrum(fn_pt);
  int N = _molecules.size();
  for(int i=0;i<N;i++)
    {
      // cout << i << " : " << endl;
      if(param_condition(_parameters[i]))
	{
	treelet_spectrum ** mol_distrib =_molecules[i]->getLabeledTreeletSpectrum();
	  for(int j=0;j<SIZE_SPECTRUM;j++)
	    {
	      for(treelet_spectrum::iterator it_mol = mol_distrib[j]->begin();
		  it_mol != mol_distrib[j]->end(); it_mol++)
		{
		  // assert((it_mol->first).length() > 0);
		  // if(j==0)
		  //   cout << MoleculeGraph::translateTreeletCode(it_mol->first) << " = "  << it_mol->second << endl;
		  if(distribution[j].find(it_mol->first) != distribution[j].end())
		      distribution[j].find(it_mol->first)->second += it_mol->second; 
		  else
		    {
		      distribution[j].insert(pair<string,double>(it_mol->first,it_mol->second)); 
		      nbTreelets++;
		    }		  
		}
	    }
	}
    }
  return distribution;
}

treelet_spectrum * MoleculesDataset::getCycleTreeletDistribution(bool (*param_condition)(double i)){
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  nbCyclesTreelets = 0;
  treelet_spectrum * distribution = new treelet_spectrum[SIZE_SPECTRUM];
  for(int i=0;i<SIZE_SPECTRUM;i++)
    distribution[i] = treelet_spectrum(fn_pt);
  int N = _molecules.size();
  for(int i=0;i<N;i++)
    {
      if(param_condition(_parameters[i]))
	{
	treelet_spectrum ** mol_distrib =_molecules[i]->getCycleTreeletSpectrum();
	  for(int j=0;j<SIZE_SPECTRUM;j++)
	    {
	      for(treelet_spectrum::iterator it_mol = mol_distrib[j]->begin();
		  it_mol != mol_distrib[j]->end(); it_mol++)
		{
		  if(distribution[j].find(it_mol->first) != distribution[j].end())
		    distribution[j].find(it_mol->first)->second += it_mol->second; 
		  else
		    {
		      distribution[j].insert(pair<string,double>(it_mol->first,it_mol->second)); 
		      nbCyclesTreelets ++;
		    }		  
		}
	    }
	}
    }
  return distribution;
}

treelet_spectrum * MoleculesDataset::getAugmentedCycleTreeletDistribution(bool (*param_condition)(double i)){
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  nbCyclesTreelets = 0;
  treelet_spectrum * distribution = new treelet_spectrum[SIZE_SPECTRUM];
  for(int i=0;i<SIZE_SPECTRUM;i++)
    distribution[i] = treelet_spectrum(fn_pt);
  int N = _molecules.size();
  for(int i=0;i<N;i++)
    {
      if(param_condition(_parameters[i]))
	{
	treelet_spectrum ** mol_distrib =_molecules[i]->getAugmentedCycleSpectrum();
	  for(int j=0;j<SIZE_SPECTRUM;j++)
	    {
	      for(treelet_spectrum::iterator it_mol = mol_distrib[j]->begin();
		  it_mol != mol_distrib[j]->end(); it_mol++)
		{
		  if(distribution[j].find(it_mol->first) != distribution[j].end())
		    distribution[j].find(it_mol->first)->second += it_mol->second; 
		  else
		    {
		      distribution[j].insert(pair<string,double>(it_mol->first,it_mol->second)); 
		      nbCyclesTreelets ++;
		    }		  
		}
	    }
	}
    }
  return distribution;
}

treelet_spectrum * MoleculesDataset::getContractedCycleTreeletDistribution(bool (*param_condition)(double i)){
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  nbContractedCyclesTreelets = 0;
  treelet_spectrum * distribution = new treelet_spectrum[SIZE_SPECTRUM];
  for(int i=0;i<SIZE_SPECTRUM;i++)
    distribution[i] = treelet_spectrum(fn_pt);
  int N = _molecules.size();
  for(int i=0;i<N;i++)
    {
      if(param_condition(_parameters[i]))
	{
	  treelet_spectrum ** mol_distrib =_molecules[i]->getContractedCycleSpectrum();
	  for(int j=0;j<SIZE_SPECTRUM;j++)
	    {
	      for(treelet_spectrum::iterator it_mol = mol_distrib[j]->begin();
		  it_mol != mol_distrib[j]->end(); it_mol++)
		{
		  if(distribution[j].find(it_mol->first) != distribution[j].end())
		    distribution[j].find(it_mol->first)->second += it_mol->second; 
		  else
		    {
		      distribution[j].insert(pair<string,double>(it_mol->first,it_mol->second)); 
		      nbContractedCyclesTreelets ++;
		    }		  
		}
	    }
	}
    }
  return distribution;
}

void MoleculesDataset::computeSpecialVectors(treelet_spectrum ** specials, int nb_treelets){
  int N = _molecules.size();
  for(int i=0;i<N;i++)
    _molecules[i]->computeSpecialLabeledTreeletVector(specials,nb_treelets);
}

cimg_library::CImg<double>* MoleculesDataset::getSpecialVector(int molecule){
  return _molecules[molecule]->getSpecialLabeledTreeletSpectrum();
}

MoleculeGraph * MoleculesDataset::getMoleculeGraph(int molecule){
  return _molecules[molecule];
}

void MoleculesDataset::setAbsParameter()
{
  for (int i=0; i<_collections.size(); ++i)
    _parameters[i]=fabs(_parameters[i]); 
}

void MoleculesDataset::setClassSign()
{
  for (int i=0; i<_collections.size(); ++i)
    {
      if(_parameters[i]<0)
	{
	  _parameters[i]=0; // classe 0 angle négatif
	}
      else
	{
	  _parameters[i]=1; // classe 1 angle positif
	}
    }
}

void MoleculesDataset::RBF(double delta)
{
  unsigned int n = _collections.size();

   for (unsigned int i=0; i<n; ++i)
	for (unsigned int j=0; j<n; ++j)
	  if (i != j)
	    _gram(i,j) = exp(- (_gram(i,i)+_gram(j,j)-2*_gram(i,j))/(2*delta*delta));
		
      for (unsigned int i=0; i<n; ++i)
	_gram (i,i)=1;
}

int MoleculesDataset::getNumberByFilename (char* name)
{
  for (unsigned int i=0;i<_filenames.size();i++)
    if(strcmp(_filenames[i],name)==0)
      return i;
  return -1;
}

double MoleculesDataset::standardDeviation(bool normalize)
{
  double max=_parameters[0];

  if(normalize)
    for (unsigned i=0; i<_parameters.size(); i++)
      if(_parameters[i]> max)
	max=_parameters[i];

  double moy=0;
  for (unsigned i=0; i<_parameters.size(); i++)
    {
      if(normalize)
	{
	  moy+=_parameters[i]/max;
	}
      else
	{
	  moy+=_parameters[i];
	}
    }
  moy=moy/_parameters.size();

  double var=0;
  for (unsigned i=0; i<_parameters.size(); i++) 
    {
      if(normalize)
	{
	   var+=(_parameters[i]/max-moy)*(_parameters[i]/max-moy);
	}
      else
	{
	  var+=(_parameters[i]-moy)*(_parameters[i]-moy);
	}
    }


  var=var/_parameters.size();

  return sqrt(var);
}



/*************************************\\
//         Version sans prise en compte de la non presence
\*************************************/

// void MoleculesDataset::computeLabeledGraphletCorrelation(){
  
//   this->computeLabeledSpectrums();
//   bool(*fn_pt)(char*,char*) = string_utils::keyComp;
//   int N = this->size();
  
//   double temp_moy = 0.0;
//   map<char *, double,  bool (*)(char *, char *)> * freq_moy = 
//     new map<char *, double,  bool (*)(char *, char *)>[SIZE_SPECTRUM];
//   for(int i=0;i<SIZE_SPECTRUM;i++)
//     freq_moy[i] = map<char *, double,  bool (*)(char *, char *)> (fn_pt);

//   //Calcul des moyennes
//   for (int i = 0; i < N ; ++i) //Parcours des molécules
//     { 
//       MoleculeGraph::spectrum_map ** spectrum = _molecules[i]->getLabeledTreeletSpectrum();
//       temp_moy += _parameters[i];
//       for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de Treelet
// 	for(MoleculeGraph::spectrum_map::iterator it = spectrum[j]->begin();
// 	    it!=spectrum[j]->end(); it ++) //Parvours de chaque Treelet
// 	  {
// 	    map<char *, double>::iterator it_freq;
// 	    if((it_freq = freq_moy[j].find(it->first))
// 	       != freq_moy[j].end())
// 	      /*Incrémentation*/
// 	      it_freq->second += it->second;
// 	    else
// 	      /*Creation*/
// 	      freq_moy[j].insert(pair<char *, double>(it->first,it->second));
// 	  }
//     }
  
//   temp_moy /= (double)N;

//   for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
//     for(map<char *, double>::iterator it_freq = freq_moy[j].begin();
// 	it_freq != freq_moy[j].end();
// 	it_freq ++)//Parcours de chaque treelet
//       it_freq->second /= (double)N;

//   //Mettre des tableaux
//   //XXX: A recoder plus proprement
//   map<char *, double,  bool (*)(char *, char *)> * correlation_num = 
//     new map<char *, double,  bool (*)(char *, char *)>[SIZE_SPECTRUM];
//   map<char *, double,  bool (*)(char *, char *)> * correlation_den_x = 
//     new map<char *, double,  bool (*)(char *, char *)>[SIZE_SPECTRUM];
//   map<char *, double,  bool (*)(char *, char *)> * correlation_den_y = 
//     new map<char *, double,  bool (*)(char *, char *)>[SIZE_SPECTRUM];
//   for(int i=0;i<SIZE_SPECTRUM;i++)
//     {
//       correlation_num[i] = map<char *, double,  bool (*)(char *, char *)> (fn_pt);
//       correlation_den_x[i] = map<char *, double,  bool (*)(char *, char *)> (fn_pt);
//       correlation_den_y[i] = map<char *, double,  bool (*)(char *, char *)> (fn_pt);
//     }

//   for (int i=0;i<N;++i)//Parcours de chaque molecule
//     { 
//       MoleculeGraph::spectrum_map ** spectrum = _molecules[i]->getLabeledTreeletSpectrum();
//       double temp = _parameters[i];

//       for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de Treelet
// 	{
// 	  for(MoleculeGraph::spectrum_map::iterator it = spectrum[j]->begin();
// 	      it!=spectrum[j]->end(); it ++) //Parvours de chaque Treelet
// 	    {
// 	      double nb_occurences = it->second;
// 	      double mean_freq = freq_moy[j].find(it->first)->second;
// 	      map<char *, double>::iterator it_tmp;
// 	      //Numerateur
// 	      if((it_tmp = correlation_num[j].find(it->first))
// 		 != correlation_num[j].end())
// 		/*Incrémentation*/
// 		it_tmp->second += (nb_occurences-mean_freq)*(temp-temp_moy);
// 	      else
// 		/*Creation*/
// 		correlation_num[j].insert(pair<char *, double>(it->first,(nb_occurences-mean_freq)*(temp-temp_moy)));

// 	      if((it_tmp = correlation_den_x[j].find(it->first))
// 		 != correlation_den_x[j].end())
// 		/*Incrémentation*/
// 		it_tmp->second += pow((nb_occurences-mean_freq),2);
// 	      else
// 		/*Creation*/
// 		correlation_den_x[j].insert(pair<char *, double>(it->first,pow((nb_occurences-mean_freq),2)));
	      
// 	      if((it_tmp = correlation_den_y[j].find(it->first))
// 		 != correlation_den_y[j].end())
// 		/*Incrémentation*/
// 		it_tmp->second += pow((temp-temp_moy),2);
// 	      else
// 		/*Creation*/
// 		correlation_den_y[j].insert(pair<char *, double>(it->first,pow((temp-temp_moy),2)));
	      
// 	    }
//       }
//     }
	
//   _labeledCorrelation = new treelet_spectrum>[SIZE_SPECTRUM];
//   for(int j=0;j<SIZE_SPECTRUM;++j)
//     {//Parcours de chaque type de Treelet
//       _labeledCorrelation[j] = treelet_spectrum> (fn_pt);
//       for(map<char*,double>::iterator it_num = correlation_num[j].begin();
// 	  it_num!=correlation_num[j].end(); it_num ++) //Parvours de chaque Treelet
// 	{
// 	  double num = (it_num->second)/N;
// 	  double den_x = correlation_den_x[j].find(it_num->first)->second;
// 	  double den_y = correlation_den_y[j].find(it_num->first)->second;
// 	  _labeledCorrelation[j].insert(pair<char *, double>(it_num->first,num/(sqrt(den_x/N)*sqrt(den_y/N))));
// 	}
//     }  
// }

/*
 * @file moleculesCorrelation.cpp
 *
 * This program computes the correlation between the Graphlet Distribution and the property og the Molecule.
 *
 * Usage : moleculesCorrelation trainset_path dataset_file
 * options:
 *
 * @author Benoit GAUZERE 
 *
 * @version 1.0.0 (2010-11-2)
 */

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "MoleculeGraph.h"
#include <pandore.h>
#include "string_utils.h"
#include "GraphletCountKernel.h"

#include "MoleculesDataset.h"

using namespace std;
using namespace pandore;

#define LIMIT_LOW_FREQUENCY 20

struct Options
{
	char* trainset_path;
	char* dataset_file;

};

/**
 * Show the program's usage.
 */

void showUsage ()
{
	cerr << "Usage : moleculesCorrelation trainset_path dataset_file" << endl;

}

/**
 * Read the command line's options.
 * 
 * @param argc The number of arguments.
 * @param argv The command line's arguments.
 * @param options The structure in which the parameters will be stored.
 */

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
	
	// The default values are fixed
	
	options->trainset_path = argv[1];
	options->dataset_file = argv[2];

}


bool comp (pair<string, double>  l, pair<string, double>  r)
{ 
  return (fabs(l.second) > fabs(r.second)); 
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
	
	/**** Unlabeled Treelet Distribution ****/
	dataset.computeGraphletCorrelation();
	double * unlabeled_distribution = dataset.getUnlabeledTreeletDistribution(param_condition);
	int nb_unlabeled_treelets = 0;
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque Graphlet
	  {
	  cout << "Graphlet " << j << " : " << dataset.getCorrelationCoeff(j) << 
	    "(" << unlabeled_distribution[j] << ")" << endl;
	  nb_unlabeled_treelets += unlabeled_distribution[j];
	  }
	// cout << "Nb Unlabeled Treelets : " << nb_unlabeled_treelets << endl;
	
	/**** Labeled Treelet Distribution ****/
	dataset.computeLabeledGraphletCorrelation();
	map<string, double,bool (*)(string, string)>* correlation = dataset.getLabeledCorrelation();
	int nb_types_treelets = 0.0;
	vector< pair<string, double> > * sort_correlations = new vector< pair<string, double> >[SIZE_SPECTRUM];
	
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  {
	    map<string, double,bool (*)(string, string)>::iterator it = correlation[j].begin();
	    for(;it != correlation[j].end();it ++)
	      {
		nb_types_treelets ++;
		sort_correlations[j].push_back(*it);
	      }	    
	  }
	
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  sort(sort_correlations[j].begin(), sort_correlations[j].end(), comp);
	
	int nb_treelets_frequents = 0;
	treelet_spectrum * distribution = dataset.getTreeletDistribution(param_condition);
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  {
	    treelet_spectrum::iterator it = distribution[j].begin();
	    for(;it != distribution[j].end();it ++)
	      {
		if(it->second > LIMIT_LOW_FREQUENCY)
		  nb_treelets_frequents ++;
	      }
	  }
	    
	double nb_treelets = 0.0;
	for(int j=0;j<SIZE_SPECTRUM;++j)//Parcours de chaque type de treelet
	  {
	    cout << "Graphlet " << j << endl; 
	    for(unsigned int i=0;i<sort_correlations[j].size();i++)
	      {
		cout << sort_correlations[j][i].first << " : " << sort_correlations[j][i].second
		     << " (" << distribution[j].find(sort_correlations[j][i].first)->second << ")"<< endl;
		nb_treelets += distribution[j].find(sort_correlations[j][i].first)->second;
	      }
	  }
	cout << "Nb Molécules : " << dataset.size() << endl;
	cout << "Nb Treelets différents : " << nb_types_treelets << endl;	
	cout << "Nb Motifs frequents : " << nb_treelets_frequents << endl;	
	cout << "Nb Treelets : " << nb_treelets << endl;
	dataset.destroy();
	
	return 0;
}

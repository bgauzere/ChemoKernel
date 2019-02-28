/*
 * @file KPrototypeEditDistance.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb  1 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cassert>
#include <algorithm>
#include <vector>
#include <float.h>
#include "KPrototypeEditDistance.h"

using namespace cimg_library;
using namespace pandore;
using namespace std;

void KPrototypeEditDistance::selectRandomPrototypes(int nb_prototypes){
  assert(nb_prototypes<=_distanceMatrix.width());
  srand (time(NULL));
  bool * flags = new bool[_distanceMatrix.width()];
  memset(flags,0,_distanceMatrix.width()*sizeof(bool));
  for(int i=0;i<nb_prototypes;i++)
    {
      int prototype = 0;
      while(flags[(prototype = rand()% _distanceMatrix.width())] == 1)
	;
      _prototypes[i] = prototype;
      flags[prototype] = 1;
    }
}

void KPrototypeEditDistance::selectCenterPrototypes(int nb_prototypes){
  CImg<double> sum_distances = CImg<double>(_distanceMatrix.width(),1,1,1,0.0);
  for(int i=0;i<_distanceMatrix.width();i++)
    for(int j=0;j<_distanceMatrix.height();j++)
      sum_distances(i) +=_distanceMatrix(i,j);
  //sum_distances(i) contains the sum of edit distances of G_i
  CImg<double> permutations = CImg<double>(_distanceMatrix.width());
  sum_distances.sort(permutations,true);
  //permutations(i) contains the initial rank of the sum, ie the index of the graph
  for(int i=0;i<nb_prototypes;i++)
    _prototypes[i] = permutations[i];
  
}

void KPrototypeEditDistance::selectMarginalPrototypes(int nb_prototypes){
  CImg<double> sum_distances = CImg<double>(_distanceMatrix.width(),1,1,1,0.0);
  for(int i=0;i<_distanceMatrix.width();i++)
    for(int j=0;j<_distanceMatrix.height();j++)
      sum_distances(i) +=_distanceMatrix(i,j);
  //sum_distances(i) contains the sum of edit distances of G_i
  CImg<double> permutations = CImg<double>(_distanceMatrix.width());
  sum_distances.sort(permutations,false);
  //permutations(i) contains the initial rank of the sum, ie the index of the graph
  for(int i=0;i<nb_prototypes;i++)
    _prototypes[i] = permutations[i];
}

void KPrototypeEditDistance::selectTargetSpherePrototypes(int nb_prototypes){
  //get Center(S)
  CImg<double> max_distances = CImg<double>(_distanceMatrix.width(),1,1,1,0.0);
  for(int i=0;i<_distanceMatrix.width();i++)
    for(int j=0;j<_distanceMatrix.height();j++)
      max_distances(i) = (_distanceMatrix(i,j) > max_distances(i))?_distanceMatrix(i,j):max_distances(i);
  
  CImg<double> permutations = CImg<double>(_distanceMatrix.width());
  max_distances.sort(permutations,true);
  int g_c = permutations[0];
  _prototypes[0] = g_c;
  //Get g_f
  double max_dist = 0.0;
  int g_f;
  for(int i=0;i<_distanceMatrix.width();i++)
    if( _distanceMatrix(g_c,i) > max_dist){
	g_f = i;
	max_dist = _distanceMatrix(g_c,i);
    }
  _prototypes[nb_prototypes-1] = g_f;;
  double w = max_dist / (double) nb_prototypes;
  bool * flags = new bool[_distanceMatrix.width()];
  memset(flags,0,_distanceMatrix.width()*sizeof(bool));
  flags[g_f] = true;
  flags[g_c] = true;
  for (int i=1;i<nb_prototypes-1;i++){
    //Find the closest to the circle of radius iw around gc
    CImg<double> distances_to_radius = CImg<double>(_distanceMatrix.width());
    for(int j=0;j<_distanceMatrix.width();j++)
      distances_to_radius(j) = fabs(_distanceMatrix(g_c,j) - i*w);
    max_distances.sort(permutations,true);
    int j=0;
    int g_i;
    while(flags[(g_i = permutations[j])] == true)
      j++;
    _prototypes[i] = g_i;
    flags[g_i] = true;
  }
}

void KPrototypeEditDistance::selectSpanningPrototypes(int nb_prototypes){
  //Get median(S)
  CImg<double> sum_distances = CImg<double>(_distanceMatrix.width(),1,1,1,0.0);
  for(int i=0;i<_distanceMatrix.width();i++)
    for(int j=0;j<_distanceMatrix.height();j++)
      sum_distances(i) +=_distanceMatrix(i,j);
  //sum_distances(i) contains the sum of edit distances of G_i
  CImg<double> permutations = CImg<double>(_distanceMatrix.width());
  sum_distances.sort(permutations,true);
  //permutations(i) contains the initial rank of the sum, ie the index of the graph
  int m = permutations[0];
  _prototypes[0] = m;
  bool * flags = new bool[_distanceMatrix.width()];
  memset(flags,0,_distanceMatrix.width()*sizeof(bool));
  flags[m] = true;
  int n = 1;
  while(_prototypes.size() < (unsigned int) nb_prototypes){
    //Get Candidates
    vector<int> candidats;
    for(int i=0;i<_distanceMatrix.width();i++)
      if(flags[i] == 0)
	candidats.push_back(i);
    
    //Compute min distances from prototypes
    CImg<double> min_distances = CImg<double>(candidats.size(),1,1,1,DBL_MAX);
    for(unsigned int c=0;c<candidats.size();++c)
      for(int i=0;i<_distanceMatrix.width();i++)
	if(flags[i] == true)
	  min_distances(c) = (_distanceMatrix(candidats[c],i) < min_distances(c))?
	    _distanceMatrix(candidats[c],i):min_distances(c);
    
    //Get farthest from prototypes
    permutations = CImg<double>(candidats.size());
    sum_distances.sort(permutations,false);
    int elu = candidats[permutations[0]];
    _prototypes[n++] = elu;
    flags[elu] = true;
  }
}

void KPrototypeEditDistance::selectKCenterPrototypes(int nb_prototypes){
  selectSpanningPrototypes(nb_prototypes);
  bool has_changed = false;
  do{
    has_changed = false;
    vector<int> * sets = new vector<int>[nb_prototypes];
    bool * flags = new bool[_distanceMatrix.width()];
    memset(flags,0,_distanceMatrix.width()*sizeof(bool));
    for(int i=0;i<nb_prototypes;i++){      
      flags[_prototypes[i]] = true;
      sets[i] = vector<int>();
      sets[i].push_back(_prototypes[i]);
    }
    for(int g=0;g<_distanceMatrix.width();g++)
      if(flags[g] == false){ //g is not a prototype
	double min_dist = DBL_MAX;
	int p_i ;
	for(int p=0;p<nb_prototypes;p++)//Find the closest prototype
	  if(_distanceMatrix(_prototypes[p],g) < min_dist){
	    min_dist = _distanceMatrix(_prototypes[p],g);
	    p_i = p;
	  }
	sets[p_i].push_back(g);
      }
    //Compute new centers
    for(int c_i=0;c_i<nb_prototypes;c_i++){
      //get Center(C_i)
      CImg<double> max_distances = CImg<double>(sets[c_i].size(),1,1,1,0.0);
      for(unsigned int i=0;i<sets[c_i].size();i++)
	for(unsigned int j=0;j<sets[c_i].size();j++){
	  double d = _distanceMatrix(sets[c_i][i],sets[c_i][j]);
	  max_distances(i) = (d > max_distances(i))?d:max_distances(i);
	}
      CImg<double> permutations = CImg<double>(sets[c_i].size());
      max_distances.sort(permutations,true);
      // /!\ Warning of many centers
      double min_max_dist = max_distances[0];
      int current_center = 0;
      bool is_prototype_stills_a_center = false;
      //Test all centers
      while(min_max_dist ==  max_distances[current_center]){
	if(sets[c_i][permutations[current_center]] == _prototypes[c_i])
	   is_prototype_stills_a_center = true;
	current_center ++;
      }
      int center = sets[c_i][permutations[0]];
      if((center != _prototypes[c_i]) && (!is_prototype_stills_a_center)){
	has_changed = true;
	_prototypes[c_i] = center;
      }
    }
    /*DEBUG*/
    for(int i=0;i<_prototypes.size();i++)
      cout << _prototypes[i] << " " ;
    cout << endl;
  }while(has_changed);
}

void KPrototypeEditDistance::selectPrototypes(int nb_prototypes, 
					      KPrototypeEditDistance::PrototypeSelectionType method){
  switch(method){
  case RandomPrototypeSelectionType: //0
    selectRandomPrototypes(nb_prototypes);
    break;
  case CenterPrototypeSelectionType: //1
    selectCenterPrototypes(nb_prototypes);
    break;
  case MarginalPrototypeSelectionType: //2
    selectMarginalPrototypes(nb_prototypes);
    break;
  case TargetSpherePrototypeSelectionType: //3
    selectTargetSpherePrototypes(nb_prototypes);
    break; 
  case SpanningPrototypeSelectionType: //4
    selectSpanningPrototypes(nb_prototypes);
    break;
  case KCenterPrototypeSelectionType: //5
    selectKCenterPrototypes(nb_prototypes);
    break;
  default:
    cerr << "Selection not defined" << endl;
    break;
  }
}

void KPrototypeEditDistance::init(const MoleculesDataset& trainset, int nb_prototypes,
			     KPrototypeEditDistance::PrototypeSelectionType method){
  _trainset = trainset;
  
  //Initialisations d,_distanceMatrix
  _d = MoleculeGraphEditDistance();
  //Compute Edit Distance Matrix
  unsigned int N = trainset.size();
  _distanceMatrix =  CImg<double>(N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      _distanceMatrix(i,j) = _d(trainset[j], trainset[i]);

  //Select Prototypes
  this->_prototypes = vector<int>(nb_prototypes);
  _prototypeDistanceMatrix = CImg<double>(trainset.size(),nb_prototypes);
  selectPrototypes(nb_prototypes,method);
  //Compute _prototypeDistanceMatrix
  for(unsigned int i=0;i<N;++i){
    _index_map.insert(pair<Collection*,int>(trainset[i],i));
    for (unsigned int j=0; j<nb_prototypes; ++j)
      _prototypeDistanceMatrix(i,j) = _distanceMatrix(i,_prototypes[j]);
  }
}

KPrototypeEditDistance::KPrototypeEditDistance(const MoleculesDataset& trainset, int nb_prototypes,
					       KPrototypeEditDistance::PrototypeSelectionType method){
  init(trainset,nb_prototypes,method);
}

KPrototypeEditDistance::KPrototypeEditDistance(const MoleculesDataset& trainset, 
					       const MoleculesDataset& testset, 
					       int nb_prototypes,
					       KPrototypeEditDistance::PrototypeSelectionType method){

  init(trainset, nb_prototypes,method);
  _prototypeDistanceMatrix.resize(trainset.size()+testset.size(),nb_prototypes);
    unsigned int N = testset.size();
    for (unsigned int i=0; i<N; i++){
    _index_map.insert(pair<Collection*,int>(testset[i],i+trainset.size()));
    for (unsigned int j=0; j<_prototypes.size(); j++){
      _prototypeDistanceMatrix(i+trainset.size(),j) = _d(trainset[j], testset[i]);
    }
  }
  
}

double KPrototypeEditDistance::operator()(pandore::Collection* c1, pandore::Collection* c2){
  if(_index_map.find(c1) == _index_map.end())
    {
      CImg<double> v(_prototypes.size());
      for (unsigned int j=0; j<_prototypes.size(); ++j)
	v(j) = _d(c1,_trainset[_prototypes[j]]);
      _prototypeDistanceMatrix.append(v.transpose());
      _index_map.insert(pair<Collection*,int>(c1,_prototypeDistanceMatrix.width()-1));
    }
  int index1 = _index_map[c1];
  
  if(_index_map.find(c2) == _index_map.end())
    {
      CImg<double> v(_prototypes.size());
      for (unsigned int j=0; j<_prototypes.size(); ++j)
	v(j) = _d(c2,_trainset[_prototypes[j]]);
      _prototypeDistanceMatrix.append(v.transpose());
      _index_map.insert(pair<Collection*,int>(c2,_prototypeDistanceMatrix.width()-1));
    }
  int index2 = _index_map[c2];
  
  double scalar =0.0;
  double norm_c1 =0.0;
  double norm_c2 =0.0;
  for(unsigned int i=0;i<_prototypes.size();i++){
    double d_1=_prototypeDistanceMatrix(index1,i);
    double d_2=_prototypeDistanceMatrix(index2,i);
    scalar += d_1*d_2;
    norm_c1 += d_1*d_1;
    norm_c2 += d_2*d_2;
  }
  // if (norm_c1*norm_c2 > 0.00000001)    
    return scalar/sqrt(norm_c1*norm_c2);
  // else
  //   return 0;
}

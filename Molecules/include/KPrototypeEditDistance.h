/*
 * @file KPrototypeEditDistance.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @TODO: Destructeur
 * @version 1.0.0 (2010-10-27)
 */

#ifndef __KPROTOTYPEEDITDISTANCE_H__
#define  __KPROTOTYPEEDITDISTANCE_H__

#include <pandore.h>
#include <deque>
#include <vector>
#include <map>
#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculeGraph.h"
#include "MoleculeGraphEditDistance.h"
#include "MoleculesDataset.h"

class KPrototypeEditDistance : public GraphKernel
{
public:
  enum PrototypeSelectionType {RandomPrototypeSelectionType, CenterPrototypeSelectionType,
			       MarginalPrototypeSelectionType,TargetSpherePrototypeSelectionType,
			       SpanningPrototypeSelectionType,KCenterPrototypeSelectionType};
  KPrototypeEditDistance(const MoleculesDataset& trainset, int nb_prototypes, 
			 PrototypeSelectionType method = RandomPrototypeSelectionType);
  KPrototypeEditDistance(const MoleculesDataset& trainset, const MoleculesDataset& testset, 
			 int nb_prototypes,PrototypeSelectionType method = RandomPrototypeSelectionType);
  void selectPrototypes(int nb_prototypes, PrototypeSelectionType method);
  double operator()(pandore::Collection* c1, pandore::Collection* c2);

private:  cimg_library::CImg<double> _distanceMatrix;//NxN matrix
  cimg_library::CImg<double> _prototypeDistanceMatrix; /* One line is a vector embedded 
							  description of a molecule of the dataset*/
  MoleculesDataset _trainset;//Keep a reference to trainset
  std::vector<int> _prototypes;
  MoleculeGraphEditDistance _d;
  std::map<pandore::Collection *, int> _index_map;/*Retrieve index from collection* */
  
  void selectRandomPrototypes(int nb_prototypes);
  void selectCenterPrototypes(int nb_prototypes);
  void selectMarginalPrototypes(int nb_prototypes);
  void selectTargetSpherePrototypes(int nb_prototypes);
  void selectSpanningPrototypes(int nb_prototypes);
  void selectKCenterPrototypes(int nb_prototypes);
  void init(const MoleculesDataset& trainset, int nb_prototypes, 
	    PrototypeSelectionType method);
};

#endif // __KPROTOTYPEEDITDISTANCE_H__

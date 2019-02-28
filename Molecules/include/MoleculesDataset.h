/**
 * @file MoleculesDataset.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 * @author Benoit Gaüzere <benoit.gauzere@insa-rouen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __MOLECULES_DATASET_H__
#define __MOLECULES_DATASET_H__

#include <deque>
#include <list>
#include <vector>
#include <pandore.h>
#include "MoleculeGraph.h"
#include "Dataset.h"
#include "CImg.h"
/**
 * @brief A dataset containing molecules with their properties.
 *
 * A dataset is composed by N molecules, with the associated properties
 * (e.g the boiling point of the molecule, or its activity), as shown below :
 *
 *		Molecule1 --> Property1
 *		Molecule2 --> Property2
 *				  ...
 *		MoleculeN --> PropertyN
 *
 * In the case of the regression, each property is a real value (e.g the boiling
 * point of the molecule).
 *
 * For the classification, a property is 0 (non selective) or 1 (selective).
 */
class MoleculesDataset;

class MoleculesDataset : public Dataset
{	
  std::deque<MoleculeGraph*> _molecules;
  double * _correlation;
  treelet_spectrum * _labeledCorrelation;
  std::deque<char *> _filenames;
  int nbTreelets;//Number of differents treelet patterns
  int nbCyclesTreelets;//Number of differents treelet patterns
  int nbContractedCyclesTreelets;//Number of differents treelet patterns extracted from _ccHypergraph
  std::vector<int> prototypes;

public:
	
  /**
   * Create an empty dataset.
   */
	
  MoleculesDataset() {}
  /**
   * Create a dataset from a given file. 
   *
   * @param path The pathname to the dataset file.
   * @param filename The file the load in the path.
   */
	
  MoleculesDataset (const char* filename);

  MoleculesDataset (const char* path, const char* filename);

  MoleculesDataset (const MoleculesDataset& d);

  /**
   * Load a dataset from a given file. 
   *
   * @param path The pathname to the dataset file.
   * @param filename The file the load in the path.
   */
	
  void loadDataset (const char* path, const char* filename);

  /*Adds to the dataset without recomputing Gram*/
  void simpleAdd (pandore::Collection* c, MoleculeGraph * m ,double parameter);

    /**
   * Reduces the size of dataset to N elements.
   *
   * Elements are randomly choosen. 
   *   *** Do it before computing Gram Matrix !! ***
   *
   * @param N The new size of dataset.
   */

  void reduceToN (unsigned int N);
  void eraseSome(std::list<int> to_erase);
  void eraseByFilename(std::vector<char *> to_erase);

  void computeSpectrums();
  void computeLabeledSpectrums();
  void computeCyclesSpectrums();
  
  void computeContractedCycleSpectrums();
  void computeAugmentedCycleSpectrums();

  void printSpectrums(); 
  void printLabeledSpectrums(); 
  void normalizeLabeledSpectrums();
  void destroy ();
  void computeGraphletCorrelation();
  void computeLabeledGraphletCorrelation();
  treelet_spectrum * getLabeledCorrelation();
  double getCorrelationCoeff(int graphlet);
  std::vector<MoleculesDataset*> splitBySize();
  void normalizeParams(double fact);
  char * getMoleculeFilename (unsigned int i) const { return _filenames[i]; }
  int getNumberByFilename (char* name);

  int getNbTreelets(){return nbTreelets;} /*Computed in getTreeletDistribution*/
  int getNbCyclesTreelets(){return nbCyclesTreelets;} /*Computed in getCycleTreeletDistribution*/

  void computeSpecialVectors(treelet_spectrum ** specials, int nbTreelets);
  cimg_library::CImg<double>* getSpecialVector(int molecule);
  int getNbBiConnexe(int molecule){return _molecules[molecule]->getNbBiConnexe();}
  std::vector< std::vector< std::pair <int, int> > > getTwoconnexeComponents(int molecule)
  { return _molecules[molecule]->getTwoconnexeComponents(); }
  void computeRelevantCycles();
  void computeCCHypergraph();
  void computeAugmentedCycles();
  void computeSimpleCycles(int k);
  std::vector<std::string> getSimpleCycleCodes(int molecule);
  
  treelet_spectrum ** getTreeletSpectrum(int molecule);
  treelet_spectrum ** getCycleTreeletSpectrum(int molecule);
  treelet_spectrum ** getContractedCycleSpectrum(int molecule);
  treelet_spectrum ** getAugmentedCycleSpectrum(int molecule);

  double * getUnlabeledTreeletDistribution(bool (*param_condition)(double i));
  treelet_spectrum * getTreeletDistribution(bool (*param_condition)(double i) );
  treelet_spectrum * getCycleTreeletDistribution(bool (*param_condition)(double i) );  
  treelet_spectrum * getAugmentedCycleTreeletDistribution(bool (*param_condition)(double i) );  
  treelet_spectrum * getContractedCycleTreeletDistribution(bool (*param_condition)(double i) );

  MoleculeGraph * getMoleculeGraph(int molecule);

  void setAbsParameter();
  void setClassSign();
  void RBF(double delta);

  /**
   * Get the standardDeviation of the properties.
   */
  double standardDeviation(bool normalize=false);

};

#endif // __MOLECULES_DATASET_H__

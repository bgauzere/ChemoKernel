/**
 * @file Dataset.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __DATASET_H__
#define __DATASET_H__

#include <deque>
#include <pandore.h>
#include <limits.h>
#include "CImg.h"
#include "GraphKernel.h"


/**
 * @brief A dataset contains objects with their properties.
 *
 * A dataset is composed by N objects, with their associated properties,
 * as shown below : 
 *
 *		Object1 --> Property1
 *		Object2 --> Property2
 *				  ...
 *		ObjectN --> PropertyN
 *
 * In the case of the regression, each property is a real value (e.g the boiling
 * point of a molecule).
 *
 * For the classification, a property is 0 (non selective) or 1 (selective).
 */
class Dataset;

class Dataset
{
protected:
  std::deque<pandore::Collection*> _collections; // The objects of the dataset
  std::deque<double> _parameters; // The parameters of the molecules
  cimg_library::CImg<double> _gram; // The Gram matrix associated to the dataset
	
public:
	
  /**
   * Create an empty dataset.
   */
	
  Dataset() {}
  
  Dataset(const Dataset& d);
  
  /**
   * Destroy the dataset.
   */
	
  virtual void destroy ();
	
  /**
   * Compute the Gram matrix associated to the dataset, with
   * the given kernel.
   * 
   * If the parameter 'normalize' is true, the Gram matrix will
   * also be normalized.
   *
   * @param kg The graph kernel to use.
   * @param normalize Normalize the matrix (default true).
   */
	
  void computeGramMatrix (GraphKernel* kg, bool normalize = true);
	
  /**
   * Get the Gram matrix associated to the dataset.
   *
   * If the parameter 'normalize' is true, the returned Gram matrix will
   * also be normalized.
   *
   * @param normalize Normalize the matrix (default true).
   *
   * @return The Gram matrix of the dataset.
   */
	
  cimg_library::CImg<double> getGramMatrix (bool normalize = true) const;
	
  /**
   * Show the Gram Matrix in the libSVM input format.
   */
	
  void showGramMatrix () const;

  /**
   * Show the Gram Matrix raw format.
   */
	
  void showGramMatrixRaw () const;
  void showGramMatrixMatlab (char * file) const;

  void showIdenticalLine();

  bool isGramMatrixPD();
  void regularizeGramMatrix();
  bool isSymmetric();
  /**
   * Delete the first object of the dataset. 
   *
   * The Gram matrix is automaticaly resized.
   */
	
  void delete_first ();

  /**
   * Delete the last object of the dataset. 
   *
   * The Gram matrix is automaticaly resized.
   */
	
  void delete_last ();

	


  /**
   * Add a new object at the end of the dataset. 
   *
   * The Gram matrix is automaticaly re-computed with the given kernel.
   *
   * @param m The object which must be added.
   * @param parameter The parameter associated to the object.
   * @param kg The graph kernel used to compute the Gram matrix.
   */
	
  void add (pandore::Collection* m, double parameter, GraphKernel* kg);

  
  /**
   * Get the size (i.e the number of objects) of the dataset.
   *
   * @return The size of the dataset.
   */
	
  unsigned int size () const { return _collections.size(); }
	
  /**
   * Get the i-th object of the dataset.
   *
   * @param i The index of the object.
   * @return The i-th object of the dataset.
   */
	
  pandore::Collection* getCollection (unsigned int i) const { return _collections[i]; }
	
  /**
   * Get the property of the i-th object of the dataset.
   *
   * @param i The index of the object.
   * @return The property of the i-th object.
   */
	
  double getParameter (unsigned int i) const { return _parameters[i]; }
	
  /**
   * Get the object which index in the dataset is i.
   *
   * @param i The index of the object.
   * @return The i-th object of the dataset.
   */
	
  pandore::Collection* operator[] (unsigned int i) const { return getCollection(i); };
	
  /**
   * Get the index of the object which has the Collection col;
   *
   * @param col The collection
   * @return The index of the object, -1 if not found.
   */
	
  int find (pandore::Collection* col) const;

  /*XXX : split deja utilisé (/home/bgauzere/dev/Molecules/Molecules/src/MoleculesDataset.cpp:48: error: no matching function for call to ‘MoleculesDataset::split(const char*, const char [2])’)*/
  //std::list<Dataset> split(int nb_classes); 

};

#endif // __DATASET_H__

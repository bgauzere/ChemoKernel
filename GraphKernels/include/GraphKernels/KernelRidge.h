/*
 * @file KernelRidge.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __KERNEL_RIDGE_H__
#define __KERNEL_RIDGE_H__

#include <vector>
#include <pandore.h>
#include "GraphKernel.h"
#include "Dataset.h"
#include "KRegression.h"

/**
 * @brief Implements the Kernel Ridge Regression.
 *
 * This class implements the Kernel Ridge regression method.
 */

class KernelRidge : public KRegression
{
  GraphKernel* _kg; // The graph kernel used
  Dataset* _dataset; // The training dataset 
  double _lambda;
	
public:
	
  /**
   * Initialize the Kernel Ridge regression method with the kernel,
   * the dataset and the value of the lambda paramter.
   *
   * @param kg The graph kernel
   * @param dataset The dataset to use
   * @param lambda The value of the lambda parameter
   */
	
  KernelRidge (GraphKernel* kg, Dataset* dataset, double lambda = 1.0)
    : _kg(kg), _dataset(dataset), _lambda(lambda)
  {
  }
  
  /**
   * Find the optimal lambda parameter of the regression for the given molecule.
   *
   * @param m The molecule to treat.
   * @param p The experimental value associated to this molecule.
   *
   * @return The optimal regression parameter for the molecule.
   */

  double optimalParameter (pandore::Collection* col, double p);
	
  /**
   * Set the lambda para<MoleculeGraph>meter of the regression.
   *
   * @param lambda The value of the lambda parameter.
   */
	
  void setLambda (double lambda) { _lambda = lambda; }
	
  /**
   * Returns the estimate value for the molecule m.
   *
   * @param m Molecule on which the regression must be computed.
   *
   * @return The estimate value for the molecule m.
   */
	
  double operator() (pandore::Collection* col);		    
};

#endif // __KERNEL_RIDGE_H__

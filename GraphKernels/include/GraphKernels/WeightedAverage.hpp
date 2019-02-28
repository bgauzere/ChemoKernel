/*
 * @file WeightedAverage.hpp
 *
 */

#ifndef __WEIGHTED_AVERAGE_H__
#define __WEIGHTED_AVERAGE_H__

#include <vector>
#include <pandore.h>
#include "GraphKernel.h"
#include "Dataset.h"
#include "KRegression.h"

/**
 * @brief Implements the Weighted Average of known values
 *
 */

class WeightedAverage : public KRegression
{
  GraphKernel* _kg; // The graph kernel used
  Dataset* _dataset; // The training dataset 
	
public:
	
  /**
   * Initialize the Weighted Average regression method with the kernel
   * and the dataset.
   *
   * @param kg The graph kernel
   * @param dataset The dataset to use
   */
	
  WeightedAverage (GraphKernel* kg, Dataset* dataset);
	
  /**
   * Returns the estimate value for the molecule m.
   *
   * @param m Molecule on which the regression must be computed.
   *
   * @return The estimate value for the molecule m.
   */
	
  double operator() (pandore::Collection* col);		    
};

#endif // __WEIGHTED_AVERAGE_H__

/**
 * @file KEditDistance.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __K_EDIT_DISTANCE_H__
#define __K_EDIT_DISTANCE_H__

#include "GraphKernel.h"
#include "EditDistance.h"
#include <pandore.h>

/**
 * @brief Implements the kernel based on the graph edit distance.
 *
 * This class implements the kernel (positive semi-definite) directly
 * based on the graph edit distance : K(G1,G2)=exp(-d(G1,G2)), where
 * d(G1,G2) is the edit distance between G1 and G2.
 */

class KEditDistance : public GraphKernel
{
  EditDistance* edit;	// The edit distance used
  
	double sigma; // The kernel parameter
  
public:
	/**
	 * Initialize the kernel with the edit distance to use.
	 *
	 * @param edit The edit distance which must be used.
	 * @param sigma The kernel sigma parameter.
	 */
	
	KEditDistance (EditDistance* edit, double sigma) : edit(edit), sigma(sigma) {}
	
	/**
	 * Returns the value of the kernel between the graphs c1 and c2.
	 *
	 * @param c1 The first graph.
	 * @param c2 The second graph.
	 *
	 * @return The value of the kernel between the two graphs.
	 */
	
	double operator() (pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __K_EDIT_DISTANCE_H__

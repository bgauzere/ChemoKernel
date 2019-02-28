/*
 * @file Kashima.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.0.0 28/03/2010
 */

#ifndef __KASHIMA_H__
#define __KASHIMA_H__

#include "TrailKernel.h"
#include <pandore.h>

using namespace pandore;

/**
 * @brief Implements the computation of the Kashima's kernel on trails.
 *
 * This class implements the computation of the Kashima's kernel on trails.
 */

class Kashima : public TrailKernel
{
 public:
	
	/**
	 * Returns the value of the kernel between the trails h1 and h2.
	 * 
	 * @param col1 The first graph
	 * @param col2 The second graph
	 * @param h1 The first trail
	 * @param h2 The second trail
	 *
	 * @return The value of the kernel between h1 and h2.
	 */
	
  double operator() (Collection* col1, Collection* col2, const trail & h1, const trail & h2);
};

#endif // __KASHIMA_H__

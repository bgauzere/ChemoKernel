/*
 * @file KRegression.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __K_REGRESSION_H__
#define __K_REGRESSION_H__

#include <pandore.h>

/**
 * @brief Base for all kernel regression methods.
 *
 * The abstract class KRegression is the base class for all
 * the kernel regression methods.
 */

class KRegression
{
public:
	
	/**
	 * Returns the estimate value for the molecule m.
	 *
	 * @param m Molecule on which the regression must be computed.
	 *
	 * @return The estimate value for the molecule m.
	 */
	
	virtual double operator() (pandore::Collection* m) = 0;
};

#endif // __K_REGRESSION_H__

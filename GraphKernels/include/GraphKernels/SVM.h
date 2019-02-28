/**
 * @file SVM.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __SVM_H__
#define __SVM_H__

#include "libsvm.h"
#include "Dataset.h"
#include "GraphKernel.h"
#include "pandore.h"
#include <deque>

/**
 * @brief Interface to the LibSVM library.
 *
 * This class is an interface to the LibSVM library (http://www.csie.ntu.edu.tw/~cjlin/libsvm/).
 *
 * It enables to perform C_SVC, NU_SVC or ONE_CLASS classification, as well as
 * EPSILON_SVR and NU_SVR regression.
 */

class SVM
{
	svm_parameter param; // The LibSVM parameters
	svm_problem prob; // The LibSVM problem
	svm_model* model; // The LibSVM model
	
	std::deque<pandore::Collection*> cols;
	std::deque<int> classes;
	
	GraphKernel* kgraph; // The graph kernel used for the prediction
	
public:
	
	enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };
	
	/**
	 * Initialize LibSVM with a trainset.
	 *
	 * @param trainset The training set.
	 * @param type The SVM type (C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR or NU_SVR).
	 * @param c The SVM parameter for C_SVC, EPSILON_SVR, and NU_SVR (default 1.0).
	 * @param nu The SVM parameter for NU_SVC, ONE_CLASS, and NU_SVR (default 0.5).
	 * @param epsilon Size of the epsilon-tube for EPSILON_SVR (default 0.1).
	 */
	
	SVM (const Dataset & trainset, GraphKernel* kgraph, int type, double c=1.0, double nu=0.5,double epsilon=0.1);
	
	/**
	 * Destroy the SVM problem.
	 */
	
	~SVM ();
	
	/**
	 * Perform a classification or a regression on the test sample.
	 *
	 * @param col The test sample.
	 *
	 * @return The predicted value for the test sample.
	 */
	
	double predict (pandore::Collection* col) const;

};

#endif // __SVM_H__

/**
 * @file KernelPCA.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __KERNEL_PCA_H__
#define __KERNEL_PCA_H__

#include <deque>
#include "pandore.h"
#include "CImg.h"
#include "GraphKernel.h"
#include "Dataset.h"

/**
 * @brief A Kernel PCA implementation.
 *
 * This class is a Kernel PCA implementation. It enables to project
 * test data onto the k-th principal eigenvector of the centered Gram matrix.
 *
 * It also enables to compute the reconstruction error in the feature space of
 * a test sample, which can be used as a measure of novelty.
 *
 * See: H. Hoffmann : Kernel PCA for Novelty Detection. In Pattern Recognition, Vol. 40, 2007, pp. 863-874
 */

class KernelPCA
{
	GraphKernel* kgraph; // The graph kernel used for kernel PCA
	std::deque<pandore::Collection*> cols; // The training set
	cimg_library::CImg<double> eigvals; // The eigenvalues of the centered Gram matrix
	cimg_library::CImg<double> eigvects; // The eigenvectors of the centered Gram matrix
	cimg_library::CImg<double> k; // The Gram matrix
	cimg_library::CImg<double> kCentered; // The centered Gram matrix
	
	double kSum;
	
public:
	
	/**
	 * Initialize the Kernel PCA method with the graph kernel
	 * that will be used.
	 *
	 * @param kgraph The graph kernel used.
	 */
	
	KernelPCA (GraphKernel* kgraph) : kgraph(kgraph) {}
	
	/**
	 * Create a Kernel PCA.
	 *
	 * @param kgraph The graph kernel used.
	 * @param dataset The trainset onto which the test sample will be projected.
	 */
	
	KernelPCA (GraphKernel* kgraph, const Dataset & dataset);
	
	/**
	 * Add a Pandore collection in the training set. If compute is set
	 * to false, the Gram matrix will not been recomputed.
	 *
	 * @param col The collection which will be added.
	 * @param compute Tell whether the Gram matrix must be computed (default false).
	 *
	 */
	
	void addToTrainingSet (pandore::Collection* col, bool compute = false);
	
	/**
	 * Compute the Gram matrix and the centered Gram matrix.
	 */
	
	void computeGramMatrix ();
	
	/**
	 * Project the training set onto the p principal
	 * eigenvectors of the centered Gram matrix.
	 *
	 * @param p The number of eigenvectors.
	 *
	 * @return The p-dimensions coordinates of the graphs in the training set.
	 */
	
	double** projectTrainingSet (unsigned int p) const;
	
	/**
	 * Project a Pandore collection onto the p principal
	 * eigenvectors of the centered Gram matrix.
	 *
	 * @param col The collection to project.
	 * @param p The number of eigenvectors.
	 *
	 * @return The p-dimensions coordinates of the projection.
	 */
	
	double* project (pandore::Collection* col, unsigned int p);
	
	/**
	 * Returns the size of the training set.
	 *
	 * @return The size of the training set.
	 */
	
	unsigned int getTrainingSetSize () const { return cols.size(); }
	
	/**
	 * Returns the reconstruction error of the given sample.
	 * This measure can then be used for novelty detection.
	 *
	 * @param col The Pandore collection.
	 * @param p The number of eigenvectors.
	 *
	 * @return The reconstruction error.
	 */
	
	double reconstructionError (pandore::Collection* col, unsigned int p);
};

#endif // __KERNEL_PCA_H__
/*
 * @file LaplacianKernel.h
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#ifndef __LAPLACIAN_KERNEL_H__
#define __LAPLACIAN_KERNEL_H__

#include <pandore.h>
#include <deque>
#include "CImg.h"
#include "GraphKernel.h"
#include "GraphEditDistance.h"
#include "Dataset.h"

#define PRINT_MATRIX 0

/**
 * @brief Implements the graph laplacian kernel.
 *
 * This class implements the Graph Laplacian Kernel. In order to obtain
 * a positive definite kernel, a regularization operator can be applied on
 * the eigen values of the normalized Laplacian matrix. 
 *
 * The possible regularization operators are :
 * 
 *      0	--  No regularization : r(eigval) = eigval
 *		1	--	Regularized Laplacian : r(eigval) = 1 + eigval*lambda^2
 *		2	--	Diffusion Process : r(eigval) = exp(1/2*eigval*lambda^2)
 *		3	--	Inverse Cosine : r(eigval) = 1/(cos(eigval*PI/4))
 *		4	--	One-step Random Walk : r(eigval) = 1/(lambda - eigval)
 */


class LaplacianKernel : public GraphKernel
{
  GraphEditDistance* edit;
	
  double sigma; // The parameter of the kernel
	
  int regularization; // The regularization type
  double lambda; // The regularization parameter
	
  cimg_library::CImg<double> K; // The gram matrix
  cimg_library::CImg<double> U_K; // The unormalized gram matrix

  cimg_library::CImg<double> W; // The weights matrix
  cimg_library::CImg<double> D; // The degree matrix
  cimg_library::CImg<double> L; // The un normalized graph Laplacian matrix
	
  std::deque<pandore::Collection*> _collections;
	
  /**
   * Compute the Graph Laplacian Kernel.
   */
  void computeKernel();

  
  //procedure de verification de la matrice inverse
  void verifInversion(cimg_library::CImg<double> matrice_inv,
		      cimg_library::CImg<double> Delta_n,
		      cimg_library::CImg<double> delta_n);

  //Normalise la matrice de Gram
  void normalizeK();
  
  double computeWeight(pandore::Collection* m1, pandore::Collection* m2);
  
  //TODO : Avirer d'ici, a mettre en static
  void printMatrix(cimg_library::CImg<double> m);
  
public:

  LaplacianKernel (const LaplacianKernel* clone);
  
	
  /**
   * Initialize the graph laplacian kernel with a trainset, the
   * testset and the value of the parameter sigma.
   *
   * @param edit The graph edit distance that will be used.
   * @param trainset The train set used to construct the graph.
   * @param testset The test set used to construct the graph.
   * @param sigma The value of the kernel parameter (default 1.0).
   * @param regularization The regularization type (default 0 : no regularization).
   * @param lambda The value of the regularization parameter (default 1.0).
   */
  LaplacianKernel (GraphEditDistance* edit, const Dataset & trainset, const Dataset & testset, double sigma = 1.0, int regularization=0, double lambda=1.0);
	
  /**
   * Initialize the graph laplacian kernel with a dataset and
   * the value of the parameter sigma.
   *
   * @param edit The graph edit distance that will be used.
   * @param trainset The train set used to construct the graph.
   * @param sigma The value of the kernel parameter (default 1.0).
   * @param regularization The regularization type (default 0 : no regularization).
   * @param lambda The value of the regularization parameter (default 1.0).
   */
	
  LaplacianKernel (GraphEditDistance* edit, const Dataset & trainset, double sigma = 1.0, int regularization=0, double lambda=1.0);

  /**
   * Create an empty graph laplacian kernel.
   *
   * @param edit The graph edit distance that will be used.
   * @param sigma The value of the kernel parameter (default 1.0).
   * @param regularization The regularization type (default 0 : no regularization).
   * @param lambda The value of the regularization parameter (default 1.0).
   */
	
  LaplacianKernel (GraphEditDistance* edit, double sigma = 1.0, int regularization=0, double lambda=1.0)
    : edit(edit), sigma(sigma), regularization(regularization), lambda(lambda) {};
	
  /**
   * Apply the graph laplacian kernel to the graphs c1 and c2.
   *
   * @param c1 The first graph.
   * @param c2 The second graph.
   *
   * @return The value of the kernel between c1 and c2.
   */
	
  double operator()(pandore::Collection* c1, pandore::Collection* c2);

  void printWeightMatrix();
  void testFastInversion(pandore::Collection * c_new);
  void fastAdd(pandore::Collection * c_new, double epsilon = 0.00000001);

};

#endif // __LAPLACIAN_KERNEL_H__

/**
 * @file Kernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Fri May  4 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __KERNEL_H__
#define __KERNEL_H__

class Kernel
{
public:
  enum KernelType {IntersectionKernelType, GaussianKernelType, InnerProductKernelType, 
		   BinaryKernelType, RandomKernelType, CompleteGaussianKernelType, 
		   InnerGaussianKernelType, PolynomialKernelType};
  virtual double operator() (double x, double y);
};

#endif // __KERNEL_H__

/**
 * @file TreeletType.h
 *
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 * 
 * @version     0.0.1 - Mer 12 Juin 2013
 *  
 * Type of kernel between treelet.
 */

#ifndef __KERNEL_TYPE_H__
#define __KERNEL_TYPE_H__

enum KernelType {IntersectionKernelType,        //0
		 GaussianKernelType,            //1
		 InnerProductKernelType,        //2
		 BinaryKernelType,              //3
		 RandomKernelType,              //4
		 CompleteGaussianKernelType,    //5
		 InnerGaussianKernelType,       //6
		 PolynomialKernelType,          //7  
		 AdaptativeGaussianKernelType,  //8 (Not implemented)     
		 JaccardKernelType              //9
  };

#endif // __KERNEL_TYPE_H__

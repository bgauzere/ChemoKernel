/**
 * @file MoleculeSpectrumKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Mar 10 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __MOLECULESPECTRUMKERNEL_H__
#define __MOLECULESPECTRUMKERNEL_H__

#include <vector>
#include <map>
#include "MoleculeGraph.h"
#include <string>

class MoleculeSpectrumKernel
{

  //Weights of each Treelet
  std::map<std::string, double,bool (*)(std::string, std::string)> _weights_treelet[SIZE_SPECTRUM];
protected:
    //List of Very Important Treelet (VIT)
  std::vector<std::string> _vit_list[SIZE_SPECTRUM]; 
  double getWeight(std::string treelet, int treelet_type);//Structure code ?
  bool isAdmitted(std::string treelet, int treelet_type);
public:
  enum KernelType {IntersectionKernelType, GaussianKernelType, InnerProductKernelType, 
		   BinaryKernelType, RandomKernelType, CompleteGaussianKernelType, 
		   InnerGaussianKernelType, PolynomialKernelType};
  // virtual ~MoleculeSpectrumKernel();
  virtual double operator() (treelet_spectrum m1, treelet_spectrum m2, 
			     int treelet_type) = 0;
  void selectTreelets(std::vector<std::string> treelets[SIZE_SPECTRUM]);
  void weightTreelets(std::map<std::string, double,bool (*)(std::string, std::string)> treelets[SIZE_SPECTRUM]);

};

#endif // __MOLECULESPECTRUMKERNEL_H__

/*
 * @file VectorFileKernel.h
 *
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 *
 * @description : Computes kernel between files representing vector :
 * mol_id v_0 v_1 ...
 * 
 * @version 1.0.0 
 */

#ifndef __VECTOR_FILE_KERNEL_H__
#define __VECTOR_FILE_KERNEL_H__

#include <pandore.h>
#include <deque>
#include <map>
#include <vector>

#include "CImg.h"
#include "GraphKernel.h"
#include "MoleculesDataset.h"
#include "Kernel.h"


class VectorFileKernel : public GraphKernel
{
  std::map<std::string, double *> vectors;//file order
  std::vector<std::string> filenames;//dataset order
  int dimension;
  int dataset_size;
  Kernel* k;
public:
  VectorFileKernel(char * vector_file, MoleculesDataset * dataset, int dimension, Kernel * k);
  void computeGramMatrices(char * gram_file);
  //check the id value in given Collection
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __GRAPH_FILE_KERNEL_H__

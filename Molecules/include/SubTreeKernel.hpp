/**
 * @file  SubTreeKernel.hpp
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 * @version   0.0.1 - Jeu 26 Sep 2013
 * 
 *  
 *  Implement the tree pattern kernel [Mahé 2009]
 * and its extension to chiral molecules [Brown 2010]
 */

#ifndef __SUBTREEKERNEL_H__
#define __SUBTREEKERNEL_H__

#include "pandore.h"

class SubTreeKernel
{

public:

  static double computeSubTreeKernel(pandore::Graph3d& g1,std::vector<std::string> nodes1,
				     std::vector<std::string> edges1,pandore::Graph3d& g2, 
				     std::vector<std::string> nodes2,std::vector<std::string> edges2,
				     int depthMax,double lambda,bool untilFlag,bool branchFlag);

  static void noTottersTransform(pandore::Graph3d& gToTransform,std::vector<std::string> nodes,
				 std::vector<std::string> edges,pandore::Graph3d & gRes,
				 std::vector<std::string> & nodesRes,std::vector<std::string> & edgesRes,
				 std::vector<bool> & trueNode);
    
};

#endif // __SUBTREEKERNEL_H__

/*-----------------------------------------------------------------------

  File        : RandomWalkKernel.h

  Description : Random walk kernel as described by Vishwanathan

  Copyright  : Francois-Xavier Dupé - http://www.greyc.ensicaen.fr/~fdupe/
               	       
  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.

  ---------------------------------------------------------------------*/

#if !defined(RANDOMWALKKERNEL_H)
#define RANDOMWALKKERNEL_H

#include "GraphKernel.h"
#include "TrailKernel.h"

/**
 * The random walk kernel
 */
class RandomWalkKernel : public GraphKernel
{
private:
  
  /**
   * The kernel on trail
   */
  TrailKernel * ktrail;

  /**
   * The descent step
   */
  double lambda;

  /**
   * The tolerance bound
   */
  double tolerance;

  /**
   * Use the tree instead of the graph
   */
  bool useTree;
  
  /**
   * The random walk kernel (unnormalize)
   */
  double walkKernel ( pandore::Collection * col1, pandore::Collection * col2 );
  
public:

  /**
   * The constructor
   * @param ktrail     the kernel on trail
   * @param lambda    the descent step
   * @param tolerance the tolerance bound for convergence
   */
  RandomWalkKernel ( TrailKernel * ktrail, double lambda, double tolerance )
    : ktrail(ktrail),lambda(lambda),tolerance(tolerance),useTree(false)
  {
  }

  /**
   * Compute the normalized kernel value between two elements
   * @param col1 the first element
   * @param col2 the second element
   * @return the kernel value
   */
  double operator() ( pandore::Collection * col1, pandore::Collection * col2 );

  /**
   * Setter for the use tree variable
   */
  void setUseTree ( bool useTree )
  {
    this->useTree = useTree;
  }
};

#endif

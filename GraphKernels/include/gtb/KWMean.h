/*-----------------------------------------------------------------------

  File        : KWMean.h

  Description : Bag of Trails Kernel based on weighted mean

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

#if !defined(KWMEAN_H)
#define KWMEAN_H

#include "BagTrailKernel.h"
#include "TrailKernel.h"
#include <queue>

class KWMean : public BagTrailKernel
{
private:

  /**
   * The parameter of the polynomial kernel
   */
  double sigma;

  /**
   * The kernel mean
   * @param col1 the graph 1
   * @param col2 the graph 2
   * @param bag1 the bag of trails 1
   * @param bag2 the bag of trails 2
   */
  double kernelWMean ( pandore::Collection * col1, pandore::Collection * col2, std::deque<trail> & bag1, std::deque<trail> & bag2 );

public:
  
  /**
   * Default constructor
   * @param btrails  the bag of trails
   * @param kern     the kernel on trails
   * @param sigma    the parameter of the polynomial kernel
   */
  KWMean ( BagOfTrails * btrails, TrailKernel * kern, double sigma );
  
  /**
   * Compute the kernel value
   * @param col1 the first graphs
   * @param col2 the second graphs
   */
  double operator() ( pandore::Collection * col1, pandore::Collection * col2 );
};

#endif

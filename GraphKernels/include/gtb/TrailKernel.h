/*-----------------------------------------------------------------------

  File        : TrailKernel.h

  Description : General Class for Trail Kernels

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

#if !defined(TRAILKERNEL_H)
#define TRAILKERNEL_H

#include "BagOfTrails.h"

class TrailKernel
{
public:
  
  /**
   * Default destructor
   */
  virtual ~TrailKernel (void) {}

  /**
   * Compute the kernel value between two trails
   * @param col1 the graph 1
   * @param col2 the graph 2
   * @param h1   the trail from graph1
   * @param h2   the trail from graph2
   * @return the kernel value
   */
  virtual double operator() ( pandore::Collection * col1, pandore::Collection * col2, const trail & h1, const trail & h2 ) = 0;
};

#endif

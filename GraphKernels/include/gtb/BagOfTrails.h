/*-----------------------------------------------------------------------

  File        : BagOfTrails.h

  Description : Interface for the bag of trails construction classes

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

#if !defined(BAGOFTRAILS_H)
#define BAGOFTRAILS_H

#include <queue>
#include <vector>
#include <pandore.h>

namespace pandore {
  class Collection;
}

/**
 * Definition of a trail
 */
struct trail
{
  /**
   * The trail (alternative sequence of nodes and edges)
   */
  std::deque<int> t;

  /**
   * The labels of the trail
   */
  std::deque<int> l;

  /**
   * The weight of the path
   */
  double weight;

  /**
   * The total weight of the container of the path
   */
  double sumWeight;
};

class BagOfTrails
{	
protected:
	std::vector<std::deque<trail> > bags;
	
 public:

  /**
   * Compute the bag of trails of 2 graphs
   * @param col the input graph
   * @return the corresponding bag of trails
   */
  virtual std::deque<trail> operator() ( pandore::Collection * col ) = 0;
};

#endif

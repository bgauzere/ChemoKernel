/*-----------------------------------------------------------------------

  File        : EditDistance.h

  Description : Approximation of the edit distance using Munkres' algorithm
  proposed in
  [RNB07] K. Riesen, M. Neuhaus and H. Bunke. Bipartite Graph Matching for
          Computing the Edit Distance of Graphs. GbRPR 2007. pp 1-12, 2007

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

#if !defined(EDITDISTANCE_H)
#define EDITDISTANCE_H

#include "CImg.h"
#include <queue>

namespace pandore
{
  class Collection;
}

class EditDistance
{
private:

  /**
   * Point structure needed for the Munkres algorithm
   */
  struct Point
  {
    int x;
    int y;
  };
  
protected:

  /**
   * Alignement structure
   */
  struct Alignement
  {
    double cost;           // The alignement cost
    std::deque<int> alX;   // Alignement on the X axis
    std::deque<int> alY;   // Alignement on the Y axis
  };
  
  /**
   * Compute the edit distance using the Munkres' algorithm
   * @param cost the cost matrix
   * @return the approximation of the edit distance
   */
  static Alignement * applyMunkres ( const cimg_library::CImg<double> & cost );
  
public:

  /**
   * Default constructor
   */
  EditDistance () { };

  /**
   * Default Destructor
   */
  virtual ~EditDistance () { };

  /**
   * Compute the edit distance between two graphs
   * @param c1 the first graph
   * @param c2 the second graph
   * @return the edit cost between the 2 graphs
    */
   virtual double operator() ( pandore::Collection * c1, pandore::Collection * c2 ) = 0;
};

#endif

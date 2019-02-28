/*-----------------------------------------------------------------------

  File        : ShapeBagTrails.cpp

  Description : Compute Bag of Trails with heuristics for shape graph

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

#if !defined(SHAPEBAGTRAILS_H)
#define SHAPEBAGTRAILS_H

#include "BagOfTrails.h"
#include "CImg.h"

class ShapeBagTrails : public BagOfTrails
{
private:

  /**
   * The maximal size of the trail into the bag of trails
   */
  int trailSize;

  /**
   * The percent of kept trails
   */
  float percent;
  
  /**
   * Filter using covering algorithm
   */
  bool covering;

  /**
   * Toggle for tree or graph trails
   */
  bool useTree;

  /**
   * Use or not the real MST
   */
  bool strictTree;

  /**
   * The number of cuts inside the histogram
   */
  int nb_cuts;

  /**
   * Compute both sigma of the divided histograms
   * @param histo     the initial histogram
   * @param partie    the selected side to calculate sigma 
   * @param cut       the last value of the cut
   * @param debut     begin of the histo
   * @param fin       end
   */
  double computeSP ( const cimg_library::CImg<int> &histo, int partie, int cut, int debut, int fin );

  /**
   * Compute the "cut"-sepration of the histogram
   * @param histo   the initial histogram
   * @param debut   the beginning index
   * @param fin     the ending index
   * @param stock   the vector of the cuts
   * @param nbcut   the number of cuts
   */
  void separate ( const cimg_library::CImg<int> &histo, int debut, int fin, double * stock, int nbcut );
  
  /**
   * Compute trails inside graph with a controlled length
   * @param col     the graph
   * @param i       the depart node
   * @param maxsize the maximal size of a trail (number of edges)
   * @return the corresponding bag of trail
   */
  std::deque<trail> getGraphTrails ( pandore::Collection * col, int i, int maxsize );
  
  /**
   * Take only the "heaviest" trails with a control on the redundancy
   * @param bag        the bag of trail
   * @param col        the graphs
   * @param redundancy the seeking redundancy before computing the optimal set
   * @return the new bag of trail
   */
  std::deque<trail> filterTrail ( std::deque<trail> & bag, pandore::Collection * col, int redundancy );

  /**
   * Take only part the "heaviest" trails
   * @param bag the bag of trail
   * @param nbp the final number of trails
   */
  std::deque<trail> percentTrail ( std::deque<trail> & bag, unsigned int nbp );

  /**
   * Compute the weight histogram of one bag and cut it into several parts
   * @param bag   the bag of trails
   * @param nbcut the number of cuts
   * @return the labelled bag of paths
   */
  std::deque<trail> histogram ( const std::deque<trail> & bag, const int nbcut );
  
public:

  /**
   * Constructor
   * @param size        the maximal size of the trails
   * @param percent     the percent of kept trails
   * @param covering    filter trails using covering algorithm with selected equilibrium
   * @param useTree     use the MST in order to compute the trails
   * @param strictTree  use the real MST not the augmented (so no loop)
   * @param nb_cuts     the number of cuts inside the histogram
   */
  ShapeBagTrails ( int size, float percent, double covering = -1.0, bool useTree = false, bool strictTree = false,
		   int nb_cuts = 0 );

  /**
   * Compute the bag of trails associated to a graph
   * @param col the graph
   * @return the corresponding bag of trails
   */
  std::deque<trail> operator() ( pandore::Collection * col );
};

#endif

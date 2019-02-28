/*-----------------------------------------------------------------------

  File        : EditDistance.cpp

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

#include "EditDistance.h"
#include <queue>
#include <iostream>

using namespace cimg_library;
using namespace std;

EditDistance::Alignement * EditDistance::applyMunkres ( const CImg<double> & cost )
{
  double finalDist = 0.0;
  unsigned int k = 0;
  CImg<bool> star  ( cost.width(), cost.height(), 1, 1, false );
  CImg<bool> prime ( cost.width(), cost.height(), 1, 1, false );
  CImg<bool> cover_col ( cost.width(), 1, 1, 1, false );
  CImg<bool> cover_lin ( cost.height(), 1, 1, 1, false );
  CImg<double> C (cost);
  double emin = 1e99;
  Point z0;

  z0.x = z0.y = 0;

  if ( C.width() > C.height() ) k = C.height();
  else k = C.width();
  
  // x = column number
  // y = row number

  // STEP 0: initialization

  for ( int y = 0; y < C.height(); y++ )
    {
      double min = 1e99;
      
      // Final smallest value
      for ( int x = 0; x < C.width(); x++ )
	if ( C(x,y) < min ) min = C(x,y);

      // Substract the smallest value
      for ( int x = 0; x < C.width(); x++ )
	C(x,y) -= min;
    }

  cimg_forXY(C,x,y)
    if ( C(x,y) == 0 )
      {
	bool tostar = true;
	
	// Look at the row
	for ( int z = 0; z < C.height(); z++ )
	  if ( (z != y) && star(x,z) )
	    {
	      tostar = false;
	      break;
	    }

	if ( !tostar ) continue;

	// Look at the column
	for ( int z = 0; z < C.width(); z++ )
	  if ( (z != x) && star(z,y) )
	    {
	      tostar = false;
	      break;
	    }
	
	if ( !tostar ) continue;

	star(x,y) = true;
      }

  // STEP 1,2,3,4:
  bool notstop = true;
  unsigned int step = 1;
  
  while ( notstop )
    {
      switch (step)
	{
	case 1: // STEP 1: covering the columns
	  {
	    for ( int x = 0; x < C.width(); x++ )
	      {
		for ( int y = 0; y < C.height(); y++ )
		  if ( star(x,y) )
		    {
		      // Cover the column
		      cover_col(x) = true;
		      break;
		    }
	      }
	    
	    unsigned int count = 0;
	    for ( int x = 0; x < C.width(); x++ )
	      if ( cover_col(x) )
		count++;
	    
	    if ( k == count ) notstop = false; // DONE
	    else step = 2; // GOTO STEP 2
	    
	    break;
	  }
	case 2: // STEP 2:
	  {
	    // Look if C contains an uncovered zero
	    bool uncovered = false;
	    for ( int x = 0; x < C.width(); x++ )
	      {
		for ( int y = 0; y < C.height(); y++ )	
		  if (( !cover_lin(y) )&&( !cover_col(x) )&&( C(x,y) == 0 ))
		    {
		      uncovered = true;
		      prime(x,y) = true;
		      z0.x = x;
		      z0.y = y;

		      bool starZero = false;
		      for ( int z = 0; z < C.width(); z++ )
			if ( star(z,y) )
			  {
			    starZero = true;

			    cover_lin(y) = true;
			    cover_col(z) = false;
	
			    step = 2; // GOTO STEP 2
			    break;
			  }

		      if ( ! starZero )
			{
			  step = 3; // GOTO STEP 3
			  break;
			}
		    }

		if ( uncovered ) break;
	      }

	    if ( ! uncovered )
	      {
		emin = 1e99;

		cimg_forXY(C,x,y)
		  if (( !cover_lin(y) )&&( !cover_col(x) )&&( C(x,y) < emin ))
		    emin = C(x,y);
		
		step = 4; // GOTO STEP 4
	      }
	    
	    break;
	  }
	case 3:
	  {
	    deque<Point> S;
	    S.push_back(z0);
	    bool zero = false; // There is a starred zero inside z0 column
	    Point z1;
	    z1.x = z1.y = 0;

	    for ( int y = 0; y < C.height(); y++ )
	      if (( C(z0.x,y) == 0 )&& star(z0.x,y) )
		{
		  zero = true;
		  z1.x = z0.x;
		  z1.y = y;
		}
		
	    while ( zero )
	      {
		S.push_back(z1);

		// Replace z0 with the prime zero in the row of z1
		for ( int x = 0; x < C.width(); x++ )
		  if (( C(x,z1.y) == 0 )&& prime(x,z1.y) )
		    {
		      z0.x = x;
		      z0.y = z1.y;
		      break;
		    }
		
		S.push_back(z0);

		zero = false;
		for ( int y = 0; y < C.height(); y++ )
		  if (( C(z0.x,y) == 0 )&& star(z0.x,y) )
		    {
		      zero = true;
		      z1.x = z0.x;
		      z1.y = y;
		    }
	      }

	    cover_lin.fill(false);
	    prime.fill(false);

	    for ( unsigned int x = 0; x < S.size(); x++ )
	      {
		if (( x % 2 ) == 0 )
		  star(S[x].x,S[x].y) = true;
		else
		  star(S[x].x,S[x].y) = false;
	      }

	    step = 1; // GOTO STEP 1
	    break;
	  }
	case 4:
	  {
	    cimg_forXY(C,x,y)
	      {
		if ( cover_lin(y) )
		  C(x,y) += emin;

		if ( !cover_col(x) )
		  C(x,y) -= emin;
	      }
	    
	    step = 2; // GOTO STEP 2
	    break;
	  }
	default:
	  {
	    cerr << "Error not a step !" << endl;
	    exit(1);
	  }
	}
    }

  // DONE:
  Alignement * align = new Alignement ();
  
  cimg_forXY(C,x,y)
    if ( star(x,y) &&( C(x,y) == 0 ))
      {
	align->alX.push_back(x);
	align->alY.push_back(y);
	finalDist += cost(x,y);
      }

  align->cost = finalDist;
  
  return align;
}

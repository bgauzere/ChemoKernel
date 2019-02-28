/*-----------------------------------------------------------------------

  File        : RandomWalkKernel.cpp

  Description : Random walk kernel as describe by Vishwanathan

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

#include "CImg.h"
#include "RandomWalkKernel.h"
#include <queue>
#include <pandore.h>

using namespace std;
using namespace pandore;
using namespace cimg_library;

double RandomWalkKernel::walkKernel ( Collection * col1, Collection * col2 )
{
  deque<int> node1; // Nodes of col1
  deque<int> node2; // Nodes of col2
  Graph3d * g1 = 0;
  Graph3d * g2 = 0;

  if ( useTree )
    {
      g1 = col1->GETPOBJECT("tree",Graph3d);
      g2 = col2->GETPOBJECT("tree",Graph3d);
    }
  else
    {
      g1 = col1->GETPOBJECT("graph",Graph3d);
      g2 = col2->GETPOBJECT("graph",Graph3d);      
    }
  
  // Use the adjancy matrix in order to optimize computations
  CImg<int> adj1 (g1->Size(),g1->Size(),1,1,-1);
  CImg<int> adj2 (g2->Size(),g2->Size(),1,1,-1);
  unsigned int nbBrch1 = 0;
  unsigned int nbBrch2 = 0;
  
  // Get the nodes with useful index
  for ( int i = 0; i < g1->Size(); ++i )
    if ( (*g1)[i] != 0 )
      {
	node1.push_back(i);
	adj1(i,i) = 0;
	GEdge * edge = (*g1)[i]->Neighbours();
	for ( ; edge != 0; edge = edge->Next() )
	  {
	    nbBrch1++;
	    adj1(i,edge->Node()) = edge->Item();
	    adj1(edge->Node(),i) = edge->Item();
	  }
      }

  for ( int i = 0; i < g2->Size(); ++i )
    if ( (*g2)[i] != 0 )
      {
	node2.push_back(i);
	adj2(i,i) = 0;
	GEdge * edge = (*g2)[i]->Neighbours();
	for ( ; edge != 0; edge = edge->Next() )
	  {
	    nbBrch2++;
	    adj2(i,edge->Node()) = edge->Item();
	    adj2(edge->Node(),i) = edge->Item();
	  }
      }

  nbBrch1 /= 2;
  nbBrch2 /= 2;
  
  // Construct the weight matrix
  int N1 = node1.size();
  int N2 = node2.size();
  CImg<double> Wx (N1*N2,N1*N2,1,1,0.0);

  for ( int i = 0; i < N1*N2; ++i )
    for ( int j = i; j < N1*N2; ++j )
      {
	trail tr1;
	trail tr2;

	int a1 = i / N2;
	int a2 = (i+N2) % N2;

	int b1 = j / N2;
	int b2 = (j+N2) % N2;
   	  
	tr1.t.push_back(node1[a1]);
	tr1.l.push_back(1);
	tr1.weight = 1.0;
	tr1.sumWeight = (double)nbBrch1;
	if ( a1 != b1 )
	  {
	    tr1.t.push_back(-adj1(node1[a1],node1[b1]));
	    tr1.t.push_back(node1[b1]);
	  }

	tr2.t.push_back(node2[a2]);
	tr2.l.push_back(1);
	tr2.weight = 1.0;
	tr2.sumWeight = (double)nbBrch2;
	if ( a2 != b2 ) 
	  {
	    tr2.t.push_back(-adj2(node2[a2],node2[b2]));
	    tr2.t.push_back(node2[b2]);
	  }

	double val = 0.0;
	
	if ( (tr1.t.size() > 1)  && (tr2.t.size() > 1)
	     && (adj1(node1[a1],node1[b1]) > -1 )
	     && (adj2(node2[a2],node2[b2]) > -1 ))
	  {
	    val = (*ktrail)(col1,col2,tr1,tr2);
	  }

	Wx(i,j) = Wx(j,i) = val;
      }
 
  // Fixed point iteration
  CImg<double> px (1,N1*N2,1,1,1.0/((double)N1*(double)N2));
  CImg<double> qx = px;
  CImg<double> x = px;
  CImg<double> xpp = px;
  double norm = 1e99;
  
  while ( norm > tolerance )
    {
      x = xpp;
      xpp = px + lambda * Wx * x;

      norm = 0.0;
      for ( int i = 0; i < x.height(); ++i )
	norm += (x(0,i)-xpp(0,i))*(x(0,i)-xpp(0,i));
    }
  
  CImg<double> val = qx.transpose() * xpp;
  
  return val(0,0);
}

double RandomWalkKernel::operator() ( Collection * col1, Collection * col2 )
{
  return walkKernel(col1,col2);
}

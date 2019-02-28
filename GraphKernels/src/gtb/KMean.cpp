/*-----------------------------------------------------------------------

  File        : KMean.cpp

  Description : Bag of Trails Kernel based the mean value of the trail comparison

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

#include "KMean.h"
#include "pandore.h"

#include <deque>
#include <cmath>

using namespace pandore;
using namespace std;

double KMean::kernelMean ( Collection * col1, Collection * col2, deque<trail> & bag1, deque<trail> & bag2 )
{
	double sum;
	
	unsigned int P1 = bag1.size();
	unsigned int P2 = bag2.size();
	
	sum = 0.0;
	
	Collection** nodes1 = (Collection**) col1->GETPARRAY("nodes", Collection);
	char* edges1 = (char*) col1->GETARRAY("edges",Char);

	Collection** nodes2 = (Collection**) col2->GETPARRAY("nodes", Collection);
	char* edges2 = (char*) col2->GETARRAY("edges",Char);
	
	for ( unsigned int i = 0; i < P1; ++i )
	{
		unsigned int s = bag1[i].t.size();
		trail t1 = bag1[i];
		
		for ( unsigned int j = 0; j < P2; ++j )
		{
			if(bag2[j].t.size() != s)
				continue;
			
			if(nodes1[t1.t[0]]->GETVALUE("atom", Char) != nodes2[bag2[j].t[0]]->GETVALUE("atom", Char))
				continue;
			
			bool equal = true;
			
			for (unsigned int k=1; k<s; k += 2)
			{
				if(edges1[-t1.t[k]] != edges2[-bag2[j].t[k]])
				{
					equal=false;
					break;
				}
				
				if(nodes1[t1.t[k+1]]->GETVALUE("atom", Char) != nodes2[bag2[j].t[k+1]]->GETVALUE("atom", Char))
				{
					equal = false;
					break;
				}
			}
			
			if (equal)
				sum += 1;
		}
	}
	
	sum  /= (double) P1 * (double)P2;

	return sum;
}

double KMean::kernelMean ( Collection * col, deque<trail> & bag )
{
	double sum;
	
	unsigned int P = bag.size();
	
	sum = 0.0;
	
	Collection** nodes = (Collection**) col->GETPARRAY("nodes",Collection);
	char* edges = (char*) col->GETARRAY("edges",Char);
	
	for ( unsigned int i = 0; i < P; ++i )
	{
		unsigned int s = bag[i].t.size();
		trail t1 = bag[i];
		
		for ( unsigned int j = i+1; j < P; ++j )
		{
			if(bag[j].t.size() != s)
				continue;
			
			if(nodes[t1.t[0]]->GETVALUE("atom", Char) != nodes[bag[j].t[0]]->GETVALUE("atom", Char))
				continue;
			
			bool equal = true;
			
			for (unsigned int k=1; k<s; k += 2)
			{
				if(edges[-t1.t[k]] != edges[-bag[j].t[k]])
				{
					equal=false;
					break;
				}
				
				if(nodes[t1.t[k+1]]->GETVALUE("atom", Char) != nodes[bag[j].t[k+1]]->GETVALUE("atom", Char))
				{
					equal = false;
					break;
				}
			}
			
			if (equal)
				sum += 1;
		}
	}
	
	sum *= 2;
	sum += P;
	sum  /= (double) P * (double)P;
	
	return sum;
}

KMean::KMean ( BagOfTrails * btrails, TrailKernel * kern)
:BagTrailKernel(btrails, kern)
{
}

double KMean::operator() ( Collection * col1, Collection * col2 )
{	
	// Construct the bagtrails
	
	deque<trail> bag1 = (*bagtrails)(col1);
	deque<trail> bag2 = (*bagtrails)(col2);

	double kerval;
	
	if (col1 == col2)
		kerval = kernelMean(col1,bag1);
	else if (bag1.size() >= bag2.size())
		kerval = kernelMean(col1,col2,bag1,bag2);
	else
		kerval = kernelMean(col2,col1,bag2,bag1);
	
	return kerval;
}

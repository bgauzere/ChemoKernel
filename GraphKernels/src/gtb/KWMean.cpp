/*-----------------------------------------------------------------------

  File        : KWMean.cpp

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

#include "KWMean.h"
#include "CImg.h"
#include <pandore.h>

using namespace std;
using namespace pandore;
using namespace cimg_library;

double KWMean::kernelWMean ( Collection * col1, Collection * col2, deque<trail> & bag1, deque<trail> & bag2 )
{
	unsigned int P1 = bag1.size();
	unsigned int P2 = bag2.size();
	
	Collection** nodes1 = (Collection**) col1->GETPARRAY("nodes", Collection);
	char* edges1 = (char*) col1->GETARRAY("edges",Char);
	
	Collection** nodes2 = (Collection**) col2->GETPARRAY("nodes", Collection);
	char* edges2 = (char*) col2->GETARRAY("edges",Char);
	
	// Compute the mean vectors
	CImg<double> w1 ( P1, 1, 1, 1, 0.0 );
	CImg<double> w2 ( P2, 1, 1, 1, 0.0 );
	double sum = 0.0;
	double v1;
	
	for ( unsigned int i = 0; i < P1; ++i )
    {
		unsigned int s1 = bag1[i].t.size();
		trail t1 = bag1[i];
		
		for ( unsigned int k = 0; k < P1; ++k )
		{
			w1(k)=0;
			
			if(bag1[k].t.size() != s1)
				continue;
			
			if(nodes1[t1.t[0]]->GETVALUE("atom", Char) != nodes1[bag1[k].t[0]]->GETVALUE("atom", Char))
				continue;
			
			bool equal = true;
			
			for (unsigned int l=1; l<s1; l += 2)
			{
				if(edges1[-t1.t[l]] != edges1[-bag1[k].t[l]])
				{
					equal=false;
					break;
				}
				
				if(nodes1[t1.t[l+1]]->GETVALUE("atom", Char) != nodes1[bag1[k].t[l+1]]->GETVALUE("atom", Char))
				{
					equal = false;
					break;
				}
			}
			
			if (equal)
				w1(k) = 1;
		}
		
		v1 = w1.sum();
		
		for ( unsigned int j = 0; j < P2; ++j )
		{
			unsigned int s2 = bag2[j].t.size();
			trail t2 = bag2[j];
			
			if(t2.t.size() != s1)
				continue;
			
			if(nodes1[t1.t[0]]->GETVALUE("atom", Char) != nodes2[t2.t[0]]->GETVALUE("atom", Char))
				continue;
			
			bool equal = true;
			
			for (unsigned int l=1; l<s1; l += 2)
			{
				if(edges1[-t1.t[l]] != edges2[-t2.t[l]])
				{
					equal=false;
					break;
				}
				
				if(nodes1[t1.t[l+1]]->GETVALUE("atom", Char) != nodes2[t2.t[l+1]]->GETVALUE("atom", Char))
				{
					equal = false;
					break;
				}
			}
			
			if (equal)
			{
				for ( unsigned int k = 0; k < P2; ++k )
				{
					w2(k)=0;
					
					if(bag2[k].t.size() != s2)
						continue;
					
					if(nodes2[t2.t[0]]->GETVALUE("atom", Char) != nodes2[bag2[k].t[0]]->GETVALUE("atom", Char))
						continue;
					
					bool equal = true;
					
					for (unsigned int l=1; l<s2; l += 2)
					{
						if(edges2[-t2.t[l]] != edges2[-bag2[k].t[l]])
						{
							equal=false;
							break;
						}
						
						if(nodes2[t2.t[l+1]]->GETVALUE("atom", Char) != nodes2[bag2[k].t[l+1]]->GETVALUE("atom", Char))
						{
							equal = false;
							break;
						}
					}
					
					if (equal)
						w2(k) = 1;
				}
				
				sum += pow(v1 * w2.sum(),sigma);
			}
		}
		
		//sum += pow(v1 * w2.sum(),sigma) * (*kern)(col1,col2,bag1[i],bag2[j]);
    }
    
	sum  /= (double) P1 * (double) P2;
	return sum;
}

KWMean::KWMean ( BagOfTrails * btrails, TrailKernel * kern, double sigma )
: BagTrailKernel(btrails,kern),sigma(sigma)
{
}

double KWMean::operator() ( Collection * col1, Collection * col2 )
{
	// Construct the bagtrails
	
	deque<trail> bag1 = (*bagtrails)(col1);
	deque<trail> bag2 = (*bagtrails)(col2);
	
	double kerval;
	
	if (bag1.size() >= bag2.size())
		kerval = kernelWMean(col1,col2,bag1,bag2);
	else
		kerval = kernelWMean(col2,col1,bag2,bag1);
	
	return kerval;
}

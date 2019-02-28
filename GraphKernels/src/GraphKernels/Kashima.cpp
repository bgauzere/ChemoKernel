/*
 * @file Kashima.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.0.0 28/03/2010
 */

#include "TrailKernel.h"
#include <pandore.h>
#include "Kashima.h"
#include <iostream>
#include <cstring>

using namespace pandore;
using namespace std;

double Kashima::operator() (Collection* col1,Collection* col2, const trail & h1, const trail & h2)
{
	if(h1.t.size() != h2.t.size())
		return 0;
	
	Collection ** nodes1 = col1->GETPARRAY("nodes",Collection);
	char* edges1 = (char*) col1->GETARRAY("edges",Char);
	Collection ** nodes2 = col2->GETPARRAY("nodes",Collection);
	char* edges2 = (char*) col2->GETARRAY("edges",Char);
	
	if(nodes1[h1.t[0]]->GETVALUE("c_atom",Char) != nodes2[h2.t[0]]->GETVALUE("c_atom",Char))
		return 0;
	
	for(unsigned int i=1; i< h1.t.size(); i += 2)
	  {
		if(edges1[-h1.t[i]] != edges2[-h2.t[i]])
			return 0;
	
		if(nodes1[h1.t[i+1]]->GETVALUE("c_atom",Char) != nodes2[h2.t[i+1]]->GETVALUE("c_atom",Char))
			return 0;
	  }
	
	return 1;
}

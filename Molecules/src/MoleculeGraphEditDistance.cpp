/*
 * @file MoleculeGraph.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "MoleculeGraphEditDistance.h"

using namespace pandore;

int MoleculeGraphEditDistance::extraNodeCost (pandore::Collection* node, pandore::GEdge* e)
{
        int p=0;
	
	while (e != NULL)
	{
		++p;
		e = e->Next();
	}
	
	return 1+p+1;
}

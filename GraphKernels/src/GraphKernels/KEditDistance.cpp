/**
 * @file KEditDistance.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */
#include "KEditDistance.h"

using namespace pandore;

double KEditDistance::operator() (Collection* c1, Collection* c2)
{	
	double dist = (*edit) (c1, c2);
	
	return exp (-dist/(2*sigma*sigma));
}

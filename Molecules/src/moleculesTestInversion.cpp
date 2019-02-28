/*
 * @file 
 *
 * 
 *
 * Usage : 
 * 
 * @author Benoit GAUZERE <benoit.gauzere@ensicaen.fr>
 * @version 1.0 (2010-10-25)
 */

#include <iostream>
 
#include <pandore.h>
#include "GraphKernel.h"
#include "KEditDistance.h"

#include "MoleculesDataset.h"
#include "LaplacianKernel.h"
#include "MoleculeGraphEditDistance.h"
using namespace std;
using namespace pandore;

int main (int argc, char** argv)
{	
        MoleculesDataset dataset (argv[1],argv[2]);
	Collection * c = dataset.getCollection(0);
	dataset.eraseSome(list<int>(1,0));
	MoleculeGraphEditDistance * edit = new MoleculeGraphEditDistance();
	LaplacianKernel * kgraph = new LaplacianKernel (edit, dataset, 1, 1, 1);
	kgraph->testFastInversion(c);
	return 0;
}


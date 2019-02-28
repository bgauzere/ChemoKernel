/* -*- c-basic-offset: 3 -*-
 *
 * PANDORE (PANTHEON Project)
 *
 * GREYC IMAGE
 * 6 Boulevard Maréchal Juin
 * F-14050 Caen Cedex France
 *
 * This file is free software. You can use it, distribute it
 * and/or modify it. However, the entire risk to the quality
 * and performance of this program is with you.
 *
 *
 * For more information, refer to:
 * http://www.greyc.ensicaen.fr/EquipeImage/Pandore/
 */

/**
 * @author Régis Clouard - 1999-10-08
 * @author Régis Clouard - 2001-04-10 (version 3.00)
 * @author Régis Clouard - 2003-01-02 (version 4.00)
 */

#include <string.h>
#include <pandore.h>
using namespace pandore;

/**
 * @file main.cpp
 * Defines some facilities for reading and writing commands arguments.
 */

/**
 * Reads a command argument list (argc, argv).
 * It follows the syntax :: parameter* [-m mask] images_in* images_out*
 * The "masking" parameter specifies if the optional mask is used to realize
 * masking = 0: neither masking nor unmasking.
 * masking = 1: masking and unmaking
 * masking = 2: only the masking operation
 * masking = 3: only unmasking operation
 */
void pandore::ReadArgs( int argc, char* argv[], int parc, int finc, int foutc, Pobject** mask,
			Pobject* objin[], Pobject* objs[], Pobject* objout[], Pobject* objd[],
			char* parv[], const char* usage, char masking ) {
   register int i,k;
   char error[255];

   // Check the number of arguments or -p option.
   // Print PROTOTYPE : (name - number of parameters - number of inputs - number of outputsk).
   if ((argc>=2)&&(!strcmp(argv[1],"-p"))) {
      sprintf(error,"%s %d %d %d",argv[0], parc,finc,foutc);
      std::cout<<error<<std::endl;
      exit(0);
   }

   // Check the number of arguments or -h option.
   // Print USAGE.
   if ((argc<=parc) || ((argc>=2)&&(!strcmp(argv[1],"-h")))) {
      sprintf(error,usage,argv[0]);
      std::cerr<<error<<std::endl;
      exit(0);
   }

   // Read parameters (all are floating point values).
   k=1;
   for (i=0;i<parc;parv[i++]=argv[k++]) ;

   // Read the Mask if any:
   // a Mask is introduced by -m flag,
   // and is supposed to be a region map.
   *mask=NULL;
   if ((argc> k) && (!strcmp("-m",argv[k]))) {
      k++; // -m
      *mask=LoadFile(((k>=argc)||(!strcmp(argv[k],"-")))? NULL : argv[k]);
      k++; // mask file
   }

   if (finc > 0) {
      // Read all input files
      for (i= 0; i<finc; i++,k++) {
	 objin[i]=LoadFile(((k>=argc)||(!strcmp(argv[k],"-")))? NULL : argv[k]);
	 if (objin[i] == NULL) {
	    Exit(FAILURE);
	 }
      }

      // If the mask is introduced, mask all input images.
      if ((*mask) && ((masking==1) || (masking==2)))
	 for (i=0; i<finc; i++)
	    objs[i]=objin[i]->Mask(*mask);
      else
	 for (i=0; i< finc; i++)
	    objs[i]=objin[i];
   }
   
   // Initialize every output objects to NULL value.
   for(i=0; i< foutc; i++)
      objd[i]=objout[i]=NULL;
}

/**
 * Write a command argument list,
 * i.e all the output object in output file.
 * masking = 0: neither masking nor unmasking.
 * masking = 1: masking and unmaking
 * masking = 2: only the masking operation
 * masking = 3: only unmasking operation
 */
void pandore::WriteArgs( int argc, char* argv[], int parc, int finc, int foutc, Pobject** mask,
	        Pobject* objin[], Pobject* [], Pobject* objout[], Pobject* objd[], char masking ) {
   int	k,i;

   // If only one output file is empty: don't save any file.
   for (i=0;i<foutc;i++)
      if (objd[i]==NULL)
	 return;

   // Search the name among output files.
   k = (*mask)? parc+1+finc + 2 : parc+1+finc;
   
   // Unfilter output image with the mask and related input image.
   // If output object is not of the same type then the related input object
   // then do not use unmasking operation for this object.
   for (i=0;i<foutc;i++) {
      if ((*mask) && ((masking==1) || (masking==3)) && (i<=finc) && (objd[i]->Type()==objin[i]->Type()))
	 objout[i]=objd[i]->UnMask(*mask,objin[i]);
      else
	 objout[i]=objd[i];
   }
   
   for (i=0;i<foutc;i++,k++) {
      if (SaveFile(objout[i],((k >= argc)||(!strcmp(argv[k],"-")))? NULL : argv[k]) == FAILURE)
 	 Exit(FAILURE);
   }
}

/**
 * Displays a normalized error message for
 * bad input images.
 * @param objin	the array of input images.
 * @param finc	the number of input images.
 */
void pandore::PrintErrorFormat( Pobject* objin[], int finc ) {
   std::cerr << "Error: input types not supported by this operator: ";
   for (int _i=0; _i<finc; _i++) {
      if (_i>0)
	 std::cerr<< " x ";
      std::cerr<< objin[_i]->Name();
   }
   std::cerr<<std::endl;
}

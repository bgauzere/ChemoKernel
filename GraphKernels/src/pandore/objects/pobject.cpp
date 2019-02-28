 /* -*- c-basic-offset: 3 ; mode: c++-*-
 *
 * PANDORE(PANTHEON Project)
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
 *For more information, refer to:
 * http://www.greyc.ensicaen.fr/EquipeImage/Pandore/
 */

/**
 * @author Régis Clouard - 1999-10-08
 * @author Francois Angot - 1999-10-08
 * @author Alexandre Duret-Lutz - 1999-10-08
 * @author Régis Clouard - 2001-04-10 (version 3)
 * @author Régis Clouard - 2003-01-02 (version 4)
 * @author Régis Clouard - 2005-10-12 (-> fixed size for header)
 */

#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <pandore.h>
using namespace pandore;

/**
 * @file pobject.cpp
 *
 */

/**
 * Inverts MSB with LSB. It depends on the size of each
 * elements.
 * @param ptr	array that contains the elements to invert.
 * @param size	size of each element
 * @param nitems	number of element to read.
 */
static void Reverse( void *ptr, size_t size,  size_t  nitems ) {
   char *pti=(char*)ptr;
   char tmp[16];
   
   for (size_t i=0;i<nitems;i++) {
      memcpy(tmp,pti,size);
      for (size_t b=0;b<size;b++)
	 *(pti++)=tmp[size-1-b];
   }
}

/**
 * Reverses bytes order of the given value.
 * Typobj is always coded with 4 bytes.
 * @param l the input value.
 */
static Typobj Reverse( Typobj l ) {
   union{
      Typobj l;
      char c[4];
   } u,v;
   
   u.l = l;
   v.c[0]=u.c[3];
   v.c[1]=u.c[2];
   v.c[2]=u.c[1];
   v.c[3]=u.c[0];
   
   return v.l;
}

/*
 * Gets the local date.
 * Date in English format: year/month/day.
 * @param data the ouput date.
 */
void pandore::Date( Char *date ) {
   time_t tb;
   struct tm *nt=NULL;
   
   tb = time( &tb );
   nt = localtime( &tb );
   // tm_year: years since 1900.
   sprintf(date,"%.3d/%.2d/%.2d", nt->tm_year-100,nt->tm_mon+1,nt->tm_mday);
}

/*
 * Creates a new object with the given magic number.
 * @param objectype	the magic number of the object to build.
 * @return the new object.
 */
Pobject *pandore::NewObject( Typobj objectype ) {
   
   Pobject *pobject;
   switch (objectype) {
   case Po_Collection : pobject= new Collection; break;
   case Po_Img1duc : pobject= new Img1duc; break;
   case Po_Img1dsl : pobject= new Img1dsl; break;
   case Po_Img1dsf : pobject= new Img1dsf; break;
   case Po_Img2duc : pobject= new Img2duc; break;
   case Po_Img2dsl : pobject= new Img2dsl; break;
   case Po_Img2dsf : pobject= new Img2dsf; break;
   case Po_Img3duc : pobject= new Img3duc; break;
   case Po_Img3dsl : pobject= new Img3dsl; break;
   case Po_Img3dsf : pobject= new Img3dsf; break;
   case Po_Reg1d : pobject= new Reg1d; break;
   case Po_Reg2d : pobject= new Reg2d; break;
   case Po_Reg3d : pobject= new Reg3d; break;
   case Po_Graph2d : pobject= new Graph2d; break;
   case Po_Graph3d : pobject= new Graph3d; break;
   case Po_Imc2duc : pobject= new Imc2duc; break;
   case Po_Imc2dsl : pobject= new Imc2dsl; break;
   case Po_Imc2dsf : pobject= new Imc2dsf; break;
   case Po_Imc3duc : pobject= new Imc3duc; break;
   case Po_Imc3dsl : pobject= new Imc3dsl; break;
   case Po_Imc3dsf : pobject= new Imc3dsf; break;
   case Po_Imx1duc : pobject= new Imx1duc; break;
   case Po_Imx1dsl : pobject= new Imx1dsl; break;
   case Po_Imx1dsf : pobject= new Imx1dsf; break;
   case Po_Imx2duc : pobject= new Imx2duc; break;
   case Po_Imx2dsl : pobject= new Imx2dsl; break;
   case Po_Imx2dsf : pobject= new Imx2dsf; break;
   case Po_Imx3duc : pobject= new Imx3duc; break;
   case Po_Imx3dsl : pobject= new Imx3dsl; break;
   case Po_Imx3dsf : pobject= new Imx3dsf; break;
   case Po_Point1d : pobject= new Point1d; break;
   case Po_Point2d : pobject= new Point2d; break;
   case Po_Point3d : pobject= new Point3d; break;
   case Po_Dimension1d : pobject= new Dimension1d; break;
   case Po_Dimension2d : pobject= new Dimension2d; break;
   case Po_Dimension3d : pobject= new Dimension3d; break;
   default : pobject = 0;
   }

   return pobject;
}

#define LOAD 0
#define SAVE 1
/*
 * Just open a Pandore file.
 * Checks filename consistency.
 */
FILE *pandore::Pfopen( const char *filename, int mode ) {
   FILE *fd;

   if (mode == PLOAD) { // LOAD
      if (filename == NULL)
	 fd = stdin;
      else if ((fd = fopen(filename,"rb")) == NULL)
	 std::cerr << "Error: can't open file: "<< filename << std::endl;
   }else{  // SAVE
      if (filename == NULL)
	 fd = stdout;
      else if ((fd = fopen(filename,"wb")) == NULL)
	 std::cerr << "Error: can't create Pandore file: "<< filename << std::endl;
   }
   return fd;
}
 
/*
 * Reads the file header.
 * Check the file format from its magic number.
 * Determines whether bytes will be inverted (mismatch LSB/MSB).
 * In case of bad format, set the potype to object (= unknown object).
 * @return	the header struct.
 */
Po_header pandore::ReadHeader( FILE *fd, int &ver, bool &invertBytes ) {
   Po_header headfile;
   
   fread((void*)&headfile,sizeof(headfile),1,fd);
   
   if (!strncmp(headfile.magic,PO_MAGIC,sizeof(headfile.magic))) { // Current version
      ver=4;
   } else if (!strncmp(headfile.magic,"PANDORE05",sizeof(headfile.magic))) {  
      ver=4;
   } else if (!strncmp(headfile.magic,"PANDORE04",sizeof(headfile.magic))) {  
      ver=4;
   }else{   
      if (!strncmp(headfile.magic,"PANDORE30",sizeof(headfile.magic)))  // Current version
	 ver=3;
      else
	 ver=0;
      headfile.Type(object);
   }
   
   // Mismatch between data and architecture (LSB/MSB).
   invertBytes=false;
   if (headfile.Type()>255) {
      headfile.Type(Reverse(headfile.Type()));
      invertBytes=true;
   }
   return headfile;
}

/**
 * Saves the file header.
 */
void pandore::SaveHeader( FILE *fd, Typobj type ) {
   Po_header headfile;
   char *name;

   // For valgrind !!!!
   headfile.unused[0]=0;
   Date(headfile.date);

   memcpy(headfile.magic,PO_MAGIC,sizeof(headfile.magic));
   memset(headfile.ident,0,sizeof(headfile.ident));
#ifndef _WIN32
   if ((name=getenv("LOGNAME")))
      strncpy(headfile.ident,name,sizeof(headfile.ident)-1);
   else
#else
   if ((name=getenv("USERNAME")))
      strncpy(headfile.ident,name,sizeof(headfile.ident)-1);
   else
#endif
      strncpy(headfile.ident,"unknown",sizeof(headfile.ident)-1);
   headfile.Type(type);
   
   fwrite(&headfile,sizeof(headfile),1,fd);
}

/**
 * Loads an object from the standard input (filename=NULL)
 * or from the file filename.
 * @param filename	the filename from which to read the object. 
 */
Pobject *pandore::LoadFile( const char *filename ) {
   FILE *fd;
   Errc result;
   Pobject *obj;
   Po_header headfile;
   int ver;
   bool invertBytes;

   if (!(fd=pandore::Pfopen(filename,PLOAD)))
      return NULL;
   
   headfile=ReadHeader(fd,ver,invertBytes);

   obj = NewObject(headfile.Type());
   result = (obj!=0) && (obj->Load(fd,invertBytes)==SUCCESS);

   if (!result) {
      if ((filename))
	 std::cerr << "Error: bad Pandore format : "<< filename << std::endl;
      else
	 std::cerr << "Error: bad Pandore format : standard input" << std::endl;
   }
   if (filename != NULL) fclose(fd);

   return (result == SUCCESS)? obj : NULL;
}

/**
 * Loads an image from a file
 * or from the standard input (filename=NULL)
 */
Errc Pobject::LoadFile( const char *filename ) {
   FILE *fd;
   Po_header headfile;
   int ver;
   bool invertBytes;

   if (!(fd=pandore::Pfopen(filename,PLOAD)))
      return FAILURE;

   headfile=ReadHeader(fd,ver,invertBytes);

   if (headfile.Type() != Type())
      return FAILURE;
   
   Errc error=Load(fd, invertBytes); //  read its own attributes and data.
   
   if (filename != NULL)
      fclose(fd);
   
   return error;
}

/*
 * Saves object obj into a file or
 * into the standard output (if filename=NULL).
 */
Errc pandore::SaveFile( const Pobject *obj, const char *filename ) {
   return obj->SaveFile(filename);
}

/*
 * Saves an object into a file
 * or into the standard output (filename=NULL).
 */
Errc Pobject::SaveFile( const char *filename ) const {
   FILE *fd;
   Errc error;
   
   if (!(fd=pandore::Pfopen(filename,PSAVE)))
      return FAILURE;
   
   pandore::SaveHeader(fd,Type());
   if ( (error=Save(fd)) == FAILURE) {
      if ((filename))
	 std::cerr << "Error: can't save Pandore file: "<< filename << std::endl;
      else
	 std::cerr << "Error: can't save Pandore file: standard output" << std::endl;
   }
   if (filename != NULL)
      fclose(fd);
   
   return error;
}

/*
 * Loads the object from the file df.
 * Object is (attributes + data).
 */
Errc Pobject::Load( FILE *df, bool invert ) {
   setInversionMode(invert);
   // Rem : Cannot use && because operator && is
   // redefined for Errc and operands are evaluted
   // in parallel ! 
   if (LoadAttributes(df)) {
      return LoadData(df);
   }
   return FAILURE;
}

/*
 * Saves the object into the file df.
 * Object is (attributes + data).
 */
Errc Pobject::Save( FILE *df ) const {
   // Rem : Cannot use && because operator && is
   // redefined for Errc and operands are evaluted
   // in parallel ! 
    if ((SaveAttributes(df)))
	return SaveData(df);
    return FAILURE;
}

/**
 * Redefinition of fread to be hardware independant.
 * @param ptr	array to store read elements.
 * @param size	size of each element
 * @param nitems	number of element to read.
 * @param stream	the stream to read. 
 */
size_t Pobject::Fdecode( void *ptr, size_t  size,  size_t  nitems,  FILE *stream ) {
   size_t ret=fread(ptr,size,nitems,stream);
   if (size > 1 && _inversionMode)
      Reverse(ptr,size,nitems);
   return ret;
}

/**
 * Redefinition of fread to be hardware independant.
 * @param ptr	array that contains the elements to write.
 * @param size	size of each element
 * @param nitems	number of element to read.
 * @param stream	the stream to read. 
 */
size_t Pobject::Fencode( const void *ptr, size_t size,  size_t  nitems, FILE *stream ) const {
   return fwrite(ptr,size,nitems,stream);
}


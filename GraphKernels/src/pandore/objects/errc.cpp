/* -*- c-basic-offset: 3; mode:c++ -*-
 *
 * PANDORE (PANTHEON Project)
 *
 * GREYC IMAGE
 * 6 Boulevard Mar√©chal Juin
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
 * @author Alexandre Duret-Lutz. 1999-10-08
 * @author RÈgis Clouard - 2001-04-10 (version 3.00)
 * @author RÈgis Clouard - 2006-09-05 (fix buf on Visual C++ 6)
 */

/**
 * @file errc.cpp
 *
 * Error class and exit values.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pandore.h>
#ifdef _MSC_VER
#include <io.h>
#endif

/**
 *
 * @param e	
 */
void pandore::Exit (const Errc &e ) { e.Exit(); }

/**
 *
 */
void pandore::Errc::Exit( ) const {
   switch (_ret) {
   case FS_RET : pandore::Exit(_val.fs); 
   case Char_RET : pandore::Exit(_val.v_Char);
   case Uchar_RET : pandore::Exit(_val.v_Uchar);
   case Short_RET : pandore::Exit(_val.v_Short);
   case Ushort_RET : pandore::Exit(_val.v_Ushort);
   case Long_RET : pandore::Exit(_val.v_Long);
   case Ulong_RET : pandore::Exit(_val.v_Ulong);
   case Float_RET : pandore::Exit(_val.v_Float);
   case Double_RET : pandore::Exit(_val.v_Double);
   }
}

/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Errc values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( FS_t statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))) {
      if (statut == FAILURE){
	 fwrite("E",sizeof(Char),1,fp);
	 fclose(fp);
	 exit(1);
      }
      else{
	 fwrite("S",sizeof(Char),1,fp);
	 fclose(fp);
	 exit(0);
      }
   }
   exit(1);
}

/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Char values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Char statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("0",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Char),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Short values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Short statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("1",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Short),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Long values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Long statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("2",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Long),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Uchar values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Uchar statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("3",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Uchar),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Ushort values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Ushort statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("4",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Ushort),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Ulong values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Ulong statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("5",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Ulong),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Float values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Float statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("6",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Float),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}
/*
 * Writes the result of an operator execution into the HOME/.pandore file.
 * This allows to handle Double values, while the exit command
 * only allows CHAR values.
 */
void pandore::Exit( Double statut ) {
   FILE *fp;
   char	nomf[256];
   char *dir;
   
   if (!(dir=getenv(HOME))) exit(2);
   strcpy(nomf,dir);
   strcat(nomf,RESFIL);
   if ((fp=fopen(nomf,"wb"))>0){
      fwrite("7",sizeof(Char),1,fp);
      fwrite(&statut,sizeof(Double),1,fp);
      fclose(fp);
      exit(0);
   }
   exit(1);
}

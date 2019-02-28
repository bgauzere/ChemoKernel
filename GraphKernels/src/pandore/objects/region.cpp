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

/*
 * @author Régis Clouard - 1999-10-08
 * @author Francois Angot - 1999-10-08
 * @author Alexandre Duret-Lutz. 1999-10-08
 * @author Régis Clouard - 2001-04-10 (version 3.00)
 * @author Régis Clouard - 2002-12-22 (version 4.00)
 */

#include <pandore.h>
using namespace pandore;

/**
 * @file region.cpp
 */

/*
 * Create a copy.
 */
Pobject *Reg1d::Clone( ) const {
   Reg1d *tmp = new Reg1d(ncol);
   *tmp = *this;
   return tmp;
}

/*
 * Load object attributes from the file.
 */
Errc Reg1d::LoadAttributes( FILE *df ) {
   Img1d<Long>::LoadAttributes(df);
   if (Fdecode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save object attributes and object data into the file.
 */
Errc Reg1d::SaveAttributes( FILE *df ) const {
   Img1d<Long>::SaveAttributes(df);
   if (Fencode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save region data.
 */
Errc Reg1d::SaveData( FILE *df ) const {
   if (Labels()<= MAXUCHAR){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x=(Uchar)*p;
	 if (Fencode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
      }
   } else if (Labels() <= MAXUSHORT){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Ushort x=(Ushort)*p;
	 if (Fencode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
      }
   } else {
      if (Fencode((void*)Vector(),sizeof(Ulong),(int)(ncol),df)<(size_t)(ncol))
	 return FAILURE;
   }
   return SUCCESS;
}

/*
 * Load object data.
 */
Errc Reg1d::LoadData( FILE *df ) {
   if (Labels()<= MAXUCHAR){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   }else if (Labels() <= MAXUSHORT){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Ushort x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   }else{
      if (Fdecode((void*)Vector(),sizeof(Ulong),(int)(ncol),df)<(size_t)(ncol))
	 return FAILURE;
   }
   return SUCCESS;
}

/*
 * This function copies val on each pixel of the `this' data.
 */
Reg1d &Reg1d::operator=( const Ulong val ) {
   Img1d<Long>::operator=(val);
   nlabels=(val>0)? val:0;
   return *this;
}


/*
 * This function copies the ims data to the `this' data.
 */
Reg1d &Reg1d::operator=( const Reg1d &rgs ) {
   nlabels = rgs.nlabels;
   memcpy(Vector(),rgs.Vector(),ncol*sizeof(Ulong));
   return *this;
}

/*
 * Create a new region map with the only unmasked labels.
 * A mask is a region map.
 * -> nlabels is still the same.
 */
Pobject *Reg1d::Mask( const Pobject *mask ) {
   if ((!mask)||(mask->Type()!=Po_Reg1d)||(((Reg1d*)mask)->Size()!=Size())){
      std::cerr << "Warning: bad mask format... ignored" << std::endl;
      return this;
   }
   Reg1d *objd = new Reg1d(Props());
   Reg1d *m=(Reg1d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();
   
   for (int i=0; i<ncol; i++, pp++,pq++,pm++)
      *pp=(*pm==0)? 0 : *pq;

   objd->Labels(Labels());
   return objd;
}

/*
 * Create a new region map with the unmasked labels and new labels.
 * A mask is a region map.
 * -> nlabels is changed to the new label count.
 * -> Labels are recompute. -> (1 .. nlabels).
 */
Pobject *Reg1d::UnMask( const Pobject *mask, const Pobject *im ) {
   if ((!mask)||(mask->Type()!=Po_Reg1d)||(((Reg1d*)mask)->Size()!=Size())||(im->Type() != Type())||(((Reg1d*)im)->Size()!=Size())){
      std::cerr << "Warning: bad unmask format... ignored" << std::endl;
      return this;
   }
   if ((mask == NULL) || (mask->Type() != Po_Reg1d) || (im->Type() != Type())){
      return this;
   }
   Reg1d *objs = (Reg1d*)im;
   Reg1d *objd = new Reg1d(Props());
   Reg1d *m=(Reg1d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();
   Ulong *ps=objs->Vector();
   
   for (int i=0; i<ncol; i++, pp++,pq++,pm++,ps++)
      *pp=(*pm==0)? *ps : *pq;

   objd->Labels(Labels());
   return objd;
}

/*
 * Create a copy.
 */
Pobject *Reg2d::Clone( ) const {
   Reg2d *tmp = new Reg2d(nrow,ncol);
   *tmp = *this;
   return tmp;
}

/*
 * Load object attributes and object data from the file.
 */
Errc Reg2d::LoadAttributes( FILE *df ) {
   Img2d<Long>::LoadAttributes(df);
   if (Fdecode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save object attributes and object data into the file.
 */
Errc Reg2d::SaveAttributes( FILE *df ) const {
   Img2d<Long>::SaveAttributes(df);
   if (Fencode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save region data.
 */
Errc Reg2d::SaveData( FILE *df ) const {
   if (Labels()<= MAXUCHAR) {
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x=(Uchar)*p;
	 if (Fencode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
      }
   } else if (Labels() <= MAXUSHORT) {
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++) {
	 Ushort x=(Ushort)*p;
	 if (Fencode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
      }
   } else {
      if (Fencode((void*)Vector(),sizeof(Ulong),(size_t)(nrow*ncol),df)<(size_t)(nrow*ncol))
	 return FAILURE;
   }

   return SUCCESS;
}

/*
 * Load object data.
 */
Errc Reg2d::LoadData( FILE *df ) {
   if (Labels()<= MAXUCHAR){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   } else if (Labels() <= MAXUSHORT){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Ushort x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   } else {
      if (Fdecode((void*)Vector(),sizeof(Ulong),(int)(nrow*ncol),df)<(size_t)(nrow*ncol))
	 return FAILURE;
   }
   return SUCCESS;
}

/*
 * This function copies val on each pixel of the `this' data.
 */
Reg2d &Reg2d::operator=( const Ulong val ) {
   Img2d<Long>::operator=(val);
   nlabels=(val>0)? val:0;
   return *this;
}

/*
 * This function copies the ims data to the `this' data.
 */
Reg2d &Reg2d::operator=( const Reg2d &ims ) {
   nlabels = ims.nlabels;
   memcpy(Vector(),ims.Vector(),nrow*ncol*sizeof(Ulong));
   return *this;
}

/*
 * Create a new region map with the only unmasked labels.
 * A mask is a region map.
 * -> nlabels is changed to the new label count.
 * -> Labels are recompute. -> (1 .. nlabels).
 */
Pobject *Reg2d::Mask( const Pobject *mask ) {
   if ((!mask)||(mask->Type()!=Po_Reg2d)||(((Reg2d*)mask)->Size()!=Size())){
      std::cerr << "Warning: bad mask format... ignored" << std::endl;
      return this;
   }

   Reg2d *objd = new Reg2d(Props());
   Reg2d *m=(Reg2d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();

   for (int i=0; i<nrow*ncol; i++, pp++,pq++,pm++)
      *pp=(*pm==0)? 0 : *pq;

   objd->Labels(Labels());
   return objd;
}

/*
 * Create a new region map with the unmasked labels and new labels.
 * A mask is a region map.
 * -> nlabels is still the same.
 */
Pobject *Reg2d::UnMask( const Pobject *mask, const Pobject *im ) {
   if ((!mask)||(mask->Type()!=Po_Reg2d)||(((Reg2d*)mask)->Size()!=Size())||(im->Type() != Type())||(((Reg2d*)im)->Size()!=Size())){
      std::cerr << "Warning: bad unmask format... ignored" << std::endl;
      return this;
   }
   if ((mask == NULL) || (mask->Type() != Po_Reg2d) || (im->Type() != Type())){
      return this;
   }

   Reg2d *objs = (Reg2d*)im;
   Reg2d *objd = new Reg2d(Props());
   Reg2d *m=(Reg2d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();
   Ulong *ps=objs->Vector();

   for (int i=0; i<nrow*ncol; i++, pp++,pq++,pm++,ps++)
      *pp=(*pm==0)? *ps : *pq;

   objd->Labels(Labels());

   return objd;
}

/*
 * Create a copy.
 */
Pobject *Reg3d::Clone( ) const {
   Reg3d *tmp = new Reg3d(ndep,nrow,ncol);
   *tmp = *this;
   return tmp;
}

/*
 * Load object attributes and object data from the file.
 */
Errc Reg3d::LoadAttributes( FILE *df ) {
   Img3d<Long>::LoadAttributes(df);
   if (Fdecode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save object attributes and object data into the file.
 */
Errc Reg3d::SaveAttributes( FILE *df ) const {
   Img3d<Long>::SaveAttributes(df);
   if (Fencode((void*)&nlabels,sizeof(nlabels),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Save object data.
 */
Errc Reg3d::SaveData( FILE *df ) const {
   if (Labels()<= MAXUCHAR){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x=(Uchar)*p;
	 if (Fencode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
      }
   } else if (Labels() <= MAXUSHORT){
	 for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	    Ushort x=(Ushort)*p;
	    if (Fencode((void*)&x,sizeof(x),1,df)<1)
	       return FAILURE;
	 }
   } else {
      size_t s = ndep*nrow*ncol;
      // Problem : Visual C++ 2005 cannot call the fwrite function to write to a buffer that is larger than 64 MB in 
      // http://support.microsoft.com/default.aspx?scid=kb;en-us;899149
#ifdef _MSC_VER
      if (s*sizeof(Ulong) < (size_t)67076095) {
#endif
	 if (Fencode((void*)Vector(),sizeof(Ulong),s,df)<s)
	    return FAILURE;
#ifdef _MSC_VER
      } else {
	 ValueType *data = Vector();
	 s = nrow*ncol;
	 for (int z=0; z<ndep; z++) {
	    if (Fencode((void*)data,sizeof(Ulong),s,df)<s)
	       return FAILURE;
	    data += s;
	 }
      }
#endif      
   }
   return SUCCESS;
}

/*
 * Load object data.
 */
Errc Reg3d::LoadData( FILE *df ) {
   if (Labels()<= MAXUCHAR){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Uchar x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   }else if (Labels() <= MAXUSHORT){
      for (Ulong *p=Vector(); p<Vector()+VectorSize();p++){
	 Ushort x;
	 if (Fdecode((void*)&x,sizeof(x),1,df)<1)
	    return FAILURE;
	 *p=(Ulong)x;
      }
   } else {
      size_t s = ndep*nrow*ncol;
      // Problem : Visual C++ 2005 cannot call the fread function to read from a buffer that is larger than 64 MB in 
      // http://support.microsoft.com/default.aspx?scid=kb;en-us;899149
#ifdef _MSC_VER
      if (s*sizeof(Ulong) < (size_t)67076095) {
#endif
	 if (Fencode((void*)Vector(),sizeof(Ulong),s,df)<s)
	    return FAILURE;
#ifdef _MSC_VER
      } else {
	 ValueType *data = Vector();
	 s = nrow*ncol;
	 for (int z=0; z<ndep; z++) {
	    if (Fencode((void*)data,sizeof(Ulong),s,df)<s)
	       return FAILURE;
	    data += s;
	 }
      }
#endif
   }
   return SUCCESS;
}

/*
 * This function copies val on each pixel of the `this' data.
 */
Reg3d &Reg3d::operator=( const Ulong val ) {
   Img3d<Long>::operator=(val);
   nlabels=(val>0)? val:0;
   return *this;
}

/*
 * This function copies the rgs data to the `this' data.
 */
Reg3d &Reg3d::operator=( const Reg3d &rgs ) {
   nlabels = rgs.nlabels;
   memcpy(Vector(),rgs.Vector(),ndep*nrow*ncol*sizeof(Ulong));
   return *this;
}

/*
 * Create a new region map with the only unmasked labels.
 * A mask is a region map.
 * -> nlabels is changed to the new label count.
 * -> Labels are recompute. -> (1 .. nlabels).
 */
Pobject *Reg3d::Mask( const Pobject *mask ) {
   if ((!mask)||(mask->Type()!=Po_Reg3d)||(((Reg3d*)mask)->Size()!=Size())){
      std::cerr << "Warning: bad mask format... ignored" << std::endl;
      return this;
   }
   Reg3d *objd = new Reg3d(Props());
   Reg3d *m=(Reg3d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();

   for(int i=0; i< ndep*nrow*ncol;i++, pq++, pm++, pp++)
      *pp= (*pm==0)? 0 : *pq;

   objd->Labels(Labels());

   return objd;
}

/*
 * Create a new region map with the unmasked labels and new labels.
 * A mask is a region map.
 * -> nlabels is changed to the new label count.
 * -> Labels are recompute. -> (1 .. nlabels).
 */
Pobject *Reg3d::UnMask( const Pobject *mask, const Pobject *im ) {
   if ((!mask)||(mask->Type()!=Po_Reg3d)||(((Reg3d*)mask)->Size()!=Size())||(im->Type() != Type())||(((Reg3d*)im)->Size()!=Size())){
      std::cerr << "Warning: bad unmask format... ignored" << std::endl;
      return this;
   }
   if ((mask == NULL) || (mask->Type() != Po_Reg3d) || (im->Type() != Type())){
      return this;
   }
   Reg3d *objs = (Reg3d*)im;
   Reg3d *objd = new Reg3d(Props());
   Reg3d *m=(Reg3d*)mask;
   Ulong *pm=m->Vector();
   Ulong *pp=objd->Vector();
   Ulong *pq=Vector();
   Ulong *ps=objs->Vector();

   for (int i=0; i<ndep*nrow*ncol; i++, pp++,pq++,pm++,ps++)
      *pp=(*pm==0)? *ps : *pq;

   objd->Labels(Labels());

   return objd;
}

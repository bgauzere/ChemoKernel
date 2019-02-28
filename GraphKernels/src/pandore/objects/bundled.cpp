/* -*- c-basic-offset: 3; mode:c++ -*-
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
 * @author Alexandre Duret-Lutz - 1999-10-08
 * @author Régis Clouard - 2001-04-10 (Version 3.0)
 * @author Régis Clouard - 2006-02-21 (Fix bug on invert read array:float)
 * @author Régis Clouard - 2006-11-10 (Fix bug on object deletion)
 */

#include <pandore.h>
using namespace pandore;

/**
 * @file bundled.cpp
 * @brief BundledObject for collections.
 */

static BundledObject *LoadBundledType( FILE *df,  const std::string &s, Long size, bool inversionMode ) {
   BundledObject *bo;

   if (s == "Char") {
      bo = new BundledValue<Char>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Char") {
	 bo = new BundledArray<Char>(new Char[size/sizeof(Char)],size/sizeof(Char),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Uchar") {
      bo = new BundledValue<Uchar>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Uchar") {
	 bo = new BundledArray<Uchar>(new Uchar[size/sizeof(Uchar)],size/sizeof(Uchar),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Short") {
      bo = new BundledValue<Short>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Short") {
	 bo = new BundledArray<Short>(new Short[size/sizeof(Short)],size/sizeof(Short),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Ushort") {
      bo = new BundledValue<Ushort>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Ushort") {
	 bo = new BundledArray<Ushort>(new Ushort[size/sizeof(Ushort)],size/sizeof(Ushort),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Long") {
      bo = new BundledValue<Long>;
      bo->Load(df,inversionMode);
   } else 
   if (s == "int") {
      bo = new BundledValue<int>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Long") {
	 bo = new BundledArray<Long>(new Long[size/sizeof(Long)],size/sizeof(Long),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Ulong") {
      bo = new BundledValue<Ulong>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Ulong") {
	 bo = new BundledArray<Ulong>(new Ulong[size/sizeof(Ulong)],size/sizeof(Ulong),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Float") {
      bo = new BundledValue<Float>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Float") {
	 bo = new BundledArray<Float>(new Float[size/sizeof(Float)],size/sizeof(Float),true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Double") {
      bo = new BundledValue<Double>;
      bo->Load(df,inversionMode);
   } else 
      if (s == "Array:Double") {
	 bo = new BundledArray<Double>(new Double[size/sizeof(Double)],size/sizeof(Double),true);
	 bo->Load(df,inversionMode);
      } else 
      return NULL;
   return bo;
}

static BundledObject *LoadBundledPobject( FILE *df,  const std::string &s, Long size, bool inversionMode ) {
   BundledObject *bo;
   Pobject **oa;

   if (s == "Pobject:Collection") {
      bo = new BundledPobject(new Collection, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Collection") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Collection * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Collection;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img1duc") {
      bo = new BundledPobject(new Img1duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img1duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img1duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img1duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img1dsl") {
      bo = new BundledPobject(new Img1dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img1dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img1dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img1dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img1dsf") {
      bo = new BundledPobject(new Img1dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img1dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img1dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img1dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img2duc") {
      bo = new BundledPobject(new Img2duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img2duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img2duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img2duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img2dsl") {
      bo = new BundledPobject(new Img2dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img2dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img2dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img2dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img2dsf") {
      bo = new BundledPobject(new Img2dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img2dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img2dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img2dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img3duc") {
      bo = new BundledPobject(new Img3duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img3duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img3duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img3duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img3dsl") {
      bo = new BundledPobject(new Img3dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img3dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img3dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img3dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Img3dsf") {
      bo = new BundledPobject(new Img3dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Img3dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Img3dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Img3dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Reg1d") {
      bo = new BundledPobject(new Reg1d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Reg1d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Reg1d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Reg1d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Reg2d") {
      bo = new BundledPobject(new Reg2d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Reg2d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Reg2d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Reg2d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Reg3d") {
      bo = new BundledPobject(new Reg3d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Reg3d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Reg3d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Reg3d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Graph2d") {
      bo = new BundledPobject(new Graph2d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Graph2d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Graph2d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Graph2d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Graph3d") {
      bo = new BundledPobject(new Graph3d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Graph3d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Graph3d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Graph3d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc2duc") {
      bo = new BundledPobject(new Imc2duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc2duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc2duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc2duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc2dsl") {
      bo = new BundledPobject(new Imc2dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc2dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc2dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc2dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc2dsf") {
      bo = new BundledPobject(new Imc2dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc2dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc2dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc2dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc3duc") {
      bo = new BundledPobject(new Imc3duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc3duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc3duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc3duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc3dsl") {
      bo = new BundledPobject(new Imc3dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc3dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc3dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc3dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imc3dsf") {
      bo = new BundledPobject(new Imc3dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imc3dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imc3dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imc3dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx1duc") {
      bo = new BundledPobject(new Imx1duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx1duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx1duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx1duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx1dsl") {
      bo = new BundledPobject(new Imx1dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx1dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx1dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx1dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx1dsf") {
      bo = new BundledPobject(new Imx1dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx1dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx1dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx1dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx2duc") {
      bo = new BundledPobject(new Imx2duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx2duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx2duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx2duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx2dsl") {
      bo = new BundledPobject(new Imx2dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx2dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx2dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx2dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx2dsf") {
      bo = new BundledPobject(new Imx2dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx2dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx2dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx2dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx3duc") {
      bo = new BundledPobject(new Imx3duc, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx3duc") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx3duc * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx3duc;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx3dsl") {
      bo = new BundledPobject(new Imx3dsl, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx3dsl") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx3dsl * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx3dsl;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Imx3dsf") {
      bo = new BundledPobject(new Imx3dsf, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Imx3dsf") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Imx3dsf * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Imx3dsf;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Point1d") {
      bo = new BundledPobject(new Point1d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Point1d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Point1d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Point1d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Point2d") {
      bo = new BundledPobject(new Point2d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Point2d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Point2d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Point2d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Point3d") {
      bo = new BundledPobject(new Point3d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Point3d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Point3d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Point3d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Dimension1d") {
      bo = new BundledPobject(new Dimension1d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Dimension1d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Dimension1d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Dimension1d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Dimension2d") {
      bo = new BundledPobject(new Dimension2d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Dimension2d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Dimension2d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Dimension2d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
   if (s == "Pobject:Dimension3d") {
      bo = new BundledPobject(new Dimension3d, true);
      bo->Load(df,inversionMode);
   } else 
      if (s == "PArray:Dimension3d") {
	 // Use 4 bytes for pointers, even with 64-bit integers architecture
	 // for compatibility reasons.
	 oa=(Pobject**)new Dimension3d * [size];
	 for (int i = 0; i< (int)(size); ++i) {
	    oa[i] = new Dimension3d;
	 }
	 bo = new BundledPArray(oa,size,true);
	 bo->Load(df,inversionMode);
      } else 
	 return NULL;
   return bo;
}

BundledObject *pandore::LoadBundledObject( FILE *df,  const std::string &s, Long size, bool inversionMode ) {
   BundledObject *bo;
   
   if ((bo=LoadBundledType(df,s,size,inversionMode)) ||
       (bo=LoadBundledPobject(df,s,size/POINTERSIZE,inversionMode))){
      if (bo->valid())
	 return bo;
      delete bo;
   }
   return NULL; 
}

/**
 * Inverts MSB with LSB. It depends on the size of each
 * elements.
 * @param ptr	array that contains the elements to invert.
 * @param size	size of each element
 * @param nitems	number of elements to read.
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
 * Redefinition of fread to be hardware independant.
 * @param ptr	array to store read elements.
 * @param size	size of each element
 * @param nitems	number of elements to read.
 * @param stream	the stream to read. 
 */
size_t BundledObject::fdecode( void *ptr, size_t  size,  size_t  nitems,  FILE *stream ) {
   size_t ret=fread(ptr,size,nitems,stream);
   if (size > 1 && _inversionMode) {
      Reverse(ptr,size,nitems);
   }
   return ret;
}

/**
 * Redefinition of fread to be hardware independant.
 * @param ptr	array that contains the elements to write.
 * @param size	size of each element
 * @param nitems	number of elements to read.
 * @param stream	the stream to read. 
 */
size_t BundledObject::fencode( void *ptr, size_t size,  size_t  nitems,FILE *stream ) const {
   return fwrite(ptr,size,nitems,stream);
}

// BUNDLEDVALUE --------------------------------------------------

template< typename T >
std::string BundledValue<T>::Type( ) const {
   return TypeName< T >::Name(); 
}

template< typename T >
Errc BundledValue<T>::Load( FILE *df , bool invert ) {
   _inversionMode=invert;
   return (_valid = !(this->fdecode((void*)&_val,sizeof(T),1,df) < 1));
}

template< typename T >
Errc BundledValue<T>::Save( FILE *df ) const {
   return (fencode((void*)&_val,sizeof(_val),1,df) < 1)?FAILURE:SUCCESS; 
}

/**
 * Returns the number of bytes  of the value (e.g., 1 for Char, 4 for Float).
 * @return	the number of bytes.
 */
template< typename T >
Long BundledValue<T>::ByteSize( ) const {
   return sizeof(T);
}

template< typename T >
BundledObject *BundledValue<T>::ToArray( ) {
   T *z = new T[1];
   *z = _val;
   return new BundledArray<T>(z,1,true);
}

template< typename T >
BundledObject *BundledValue<T>::Clone( ) const { 
   return new BundledValue( _val );
}

// BUNDLEDARRAY --------------------------------------------------

template< typename T >
BundledArray<T>::BundledArray( T *array, Long nbelt, bool allocated ) : BundledObject(), _val(array), _s(nbelt*sizeof(T)), _allocated(allocated) { }

template< typename T >
BundledArray<T>::~BundledArray() {
   if (_allocated)
      delete[] _val;
}

template< typename T >
BundledObject *BundledArray<T>::Clone( ) const {
   T *t = new T[_s/sizeof(T)];
   for (int i = 0; i < (int)(_s/sizeof(T)); ++i)
      t[i] = _val[i];
   return new BundledArray( t, _s/sizeof(T),true );
}

template< typename T >
void BundledArray<T>::Append( BundledObject *bo ) {
   BundledArray *ba = dynamic_cast< BundledArray * >(bo);
   if (!ba)
      return;
   T *tmp = new T[_s/sizeof(T) + ba->NbrElements()];
   int i,j;
   for (i = 0; i < (int)(_s/sizeof(T)) ; ++i)
      tmp[i] = _val[i];
   for (j = 0; j < ba->NbrElements(); ++j)
      tmp[i+j] = ba->Array()[j];
   if (_allocated) {
      delete _val;
   }
   _s += ba->NbrElements()*sizeof(T);
   _val = tmp;
}
  
template< typename T >
std::string BundledArray<T>::Type( ) const {
   return "Array:" + TypeName< T >::Name(); 
}

template< typename T >
Errc BundledArray<T>::Load( FILE *df, bool invert ) {
   _inversionMode=invert;
   _valid = !(this->fdecode((void*)_val,sizeof(T),_s/sizeof(T),df) < 1);
   return _valid;
}

template< typename T >
Errc BundledArray<T>::Save( FILE *df ) const{ 
   return (fencode((void*)_val,_s/sizeof(T),sizeof(T),df) < 1)?FAILURE:SUCCESS; 
}
  
template< typename T >
Long BundledArray<T>::ByteSize( ) const {
   return _s;
}

template< typename T >
Long BundledArray<T>::NbrElements( ) const {
   return _s/sizeof(T);
}

// POBJECT  --------------------------------------------------

BundledPobject::BundledPobject( Pobject *po , bool allocated ) : BundledObject(),_val(po), _allocated(allocated) { }

BundledPobject::~BundledPobject() {
   if (_allocated)
      delete _val;
}

/*
 * Clones a Pandore object.
 */
BundledObject *BundledPobject::Clone( ) const {
   return new BundledPobject( (_val->Clone()), true);
}

std::string BundledPobject::Type( ) const {
   return "Pobject:" + _val->Name(); 
}
  
Errc BundledPobject::Load( FILE *df, bool invert ) {
   _inversionMode=invert;
   return (_valid= (_val->Load(df,invert)==SUCCESS));
}

Errc BundledPobject::Save( FILE *df ) const {
   return (_val->Save(df)); 
}

Long BundledPobject::ByteSize( ) const {
   return ((Imx3duc*)_val)->VectorSize();
}

BundledObject *BundledPobject::ToArray( ) {
   Pobject **z = new Pobject*[1];
   *z = _val;
   return new BundledPArray(z,1,true);
}

// PARRAY  --------------------------------------------------

BundledPArray::BundledPArray( Pobject **array, Long nbelt, bool allocated ): BundledObject(), _val(array), _s(nbelt*sizeof(Pobject*)), _allocated(allocated) { }

BundledPArray::~BundledPArray() {
   if (_allocated) {
      for (int i = 0 ; i < (int)(_s/sizeof(Pobject*)); ++i)
	 delete _val[i];
      delete[] _val;
   }
}

/*
 * Clones an array of Pandore Objects.
 */
BundledObject *BundledPArray::Clone( ) const {
   Pobject **t = new Pobject*[_s/sizeof(Pobject*)];
   for (int i = 0; i < (int)(_s/sizeof(Pobject*)); ++i)
      t[i] = (Pobject*)(_val[i]->Clone());
   return new BundledPArray( t, _s/sizeof(Pobject*),true);
}

void BundledPArray::Append( BundledObject *bo ) {
   BundledPArray *ba = dynamic_cast<BundledPArray*>(bo);
   if (!ba)
      return;
   Pobject **tmp = new Pobject*[_s/sizeof(Pobject*) + ba->NbrElements()];
   int i, j;
   for (i = 0; i < (int)(_s/sizeof(Pobject*)) ; ++i)
      tmp[i] = (Pobject*)(_val[i]->Clone());
   for (j = 0; j < ba->NbrElements(); ++j)
      tmp[i+j] =(Pobject*)(ba->PArray()[j]->Clone());
   if (_allocated) {
      for (int i = 0 ; i < (int)(_s/sizeof(Pobject*)); ++i)
	 delete _val[i];
      delete _val;
   }
   _val = tmp;
   _s += ba->NbrElements()*sizeof(Pobject*);
}
  
std::string BundledPArray::Type( ) const {
   return "PArray:" + _val[0]->Name(); 
}

Errc BundledPArray::Load( FILE *df, bool invert ) {
   _inversionMode=invert;
   _valid = true;
   for (int i = 0; i< (int)(_s/sizeof(Pobject*)) && _valid; ++i) {
      if (_val[i]->Load(df,invert) == FAILURE)
	 _valid = false;
   }
   return _valid;
}

Errc BundledPArray::Save( FILE *df ) const {
   for (int i = 0; i< (int)(_s/sizeof(Pobject*)); ++i)
      if (_val[i]->Save(df) == FAILURE)
	 return FAILURE;
   return SUCCESS; 
}

Long BundledPArray::ByteSize( ) const {
   return _s;
}

Long BundledPArray::NbrElements( ) const {
   return _s/sizeof(Pobject*);
}

// --------------------------------------------------

template class pandore::BundledValue<Char>;
template class pandore::BundledArray<Char>;
template class pandore::BundledValue<Uchar>;
template class pandore::BundledArray<Uchar>;
template class pandore::BundledValue<Short>;
template class pandore::BundledArray<Short>;
template class pandore::BundledValue<Ushort>;
template class pandore::BundledArray<Ushort>;
template class pandore::BundledValue<Long>;
template class pandore::BundledArray<Long>;
template class pandore::BundledValue<Ulong>;
template class pandore::BundledArray<Ulong>;
template class pandore::BundledValue<Float>;
template class pandore::BundledArray<Float>;
template class pandore::BundledValue<Double>;
template class pandore::BundledArray<Double>;


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
 * http://www.greyc.ensicaen.fr/~regis
 */

/**
 * @author Alexandre Duret-Lutz - 1999-10-08
 * @author Régis Clouard - 2001-04-10 (version 3.00)
 * @author Régis Clouard - 2006-02-20 (comptability little/big endian)
 * @author Régis Clouard - 2007-01-15 (fix bug on Delete)
 */

/**
 * @file collection.cpp
 *
 * Collection of objects
 */

#include <pandore.h>
using namespace pandore;

//void* pandore::__dummy;

void Collection::Set( const std::string &name, BundledObject *bo ) {
   std::map< std::string, BundledObject * >::iterator i = _objs.find(name);
   if (i != _objs.end()) {
      delete i->second;
      i->second = bo;
   } else	
      _objs[name] = bo;
}

BundledObject *Collection::Get( const std::string &name ) const {
   std::map< std::string, BundledObject * >::const_iterator i = _objs.find(name);
   return (i != _objs.end())?i->second:NULL;
}

void Collection::Erase( const std::string &name ) {
   std::map< std::string, BundledObject * >::iterator i = _objs.find(name);
   
   if ( i == _objs.end() )
      return;  
   delete i->second;
   _objs.erase(i);
}

/**
 * Bug: map cannot be desallocated item by item with erase(i).
 *      So, we use clear().
 */
void Collection::Delete( ) {
   std::map< std::string, BundledObject * >::iterator i;
   for (i = _objs.begin(); i != _objs.end(); ++i) {
      delete i->second;
   }
   _objs.clear();
}

bool Collection::Exists( const std::string &name ) const {
   std::map< std::string, BundledObject * >::const_iterator i = _objs.find(name);
   return i != _objs.end();
}

void Collection::Rename( const std::string &from, const std::string &to ) {
   std::map< std::string, BundledObject * >::iterator i = _objs.find(from);
   
   if ( i == _objs.end() )
      return;  
   _objs[to] = i->second;
   _objs.erase(i);
}

std::string Collection::GetType( const std::string &name ) const {
   BundledObject *bo = Get(name);
   if (!bo) { // Panic
      std::cerr << "Error: no attribute `"<< name.c_str() <<"' in collection." << std::endl;
      Exit(FAILURE);
   }
   return bo->Type();
}

Pobject *Collection::Clone( ) const {
   std::map< std::string, BundledObject * >::const_iterator i;
   Collection *tmp = new Collection;
   for (i = _objs.begin(); i != _objs.end(); ++i) {
//      std::cerr<< "Collection :: Clone object " <<i->first<<std::endl;
      tmp->Set( i->first, i->second->Clone() );
   }
   return tmp;
}

/*
 * rem: fread "size" is not in LoadAttributes
 * because size is not an attribute
 * but a computed value.
 */
Errc Collection::LoadData( FILE *file ) {
   Long size;
   Long attr[3];

   // number of elements
   if (Fdecode((void*)&size,sizeof(size),1,file) < 1)
      return FAILURE;

   for ( int i = 0; i<size; ++i ) {   
      // bytesize of the element.
      if (Fdecode((void*)attr,sizeof(attr[0]),3,file) < 1)
	 return FAILURE;

      char *type = new char[attr[1]+1];
      char *name = new char[attr[2]+1];

      // Important: keep fread versus Fdecode
      // for chars.
      if ((fread(type,attr[1],1,file) < 1) ||
	  (fread(name,attr[2],1,file) < 1))
	 return FAILURE;

      type[attr[1]] = 0;
      name[attr[2]] = 0;
      BundledObject *bo = LoadBundledObject(file,type,attr[0],_inversionMode);
      if (!bo) {
	 std::cerr <<"ERROR: problem during reading ..." << std::endl;
	 return FAILURE;
      }
      Set(name,bo);

      delete[] name;
      delete[] type;
   }

   return SUCCESS;
}

Errc Collection::SaveData( FILE *file ) const {
   std::map< std::string, BundledObject * >::const_iterator i;
   Long attr[3];
   std::string s;
   
   for (i = _objs.begin(); i!= _objs.end(); ++i) {
      s = i->second->Type();
      // case of pobject array: use 4 bytes for pointers.
      if ( s.find("PArray:") != std::string::npos ) {
	 attr[0] = (Long)(i->second->ByteSize()*POINTERSIZE)/sizeof(void*);
      } else {
	 attr[0] = (Long)i->second->ByteSize();
      }
      attr[1] = (Long)s.size();
      attr[2] = (Long)i->first.size();
      
      if ((Fencode((void*)attr,sizeof(attr),1,file) < 1) ||
	  (Fencode((void*)s.c_str(),attr[1],1,file) < 1)        ||
	  (Fencode((void*)i->first.c_str(),attr[2],1,file) < 1) ||
	  (i->second->Save(file) != SUCCESS)) {
	 return FAILURE;
      }
   }
   
   return SUCCESS;
}

std::list< std::string > Collection::List( ) const  {
   std::list< std::string > l;
   std::map< std::string, BundledObject * >::const_iterator i;
   
   for (i = _objs.begin(); i!=_objs.end(); ++i)
      l.push_back(i->first);
   return l;
}

Errc Collection::NbOf( const std::string &s, std::string &type_out, Long &number_out, Long &minsize_out ) const {
   char name[255];
   Long nbrcomp = 0;
   Long min = MAXLONG;
   
   sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);	  
   std::string type = GetType(name);
   if (type == "Array:Char") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Char);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Uchar") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Uchar);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Short") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Short);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Ushort") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Ushort);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Long") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Long);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Ulong") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Ulong);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Float") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Float);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
   if (type == "Array:Double") {
      do {
	 int n_ = (int)this->GETARRAYSIZE(name,Double);
	 if (n_ < min)
	    min = n_;
	 ++nbrcomp;
	 sprintf(name, "%s.%d", s.c_str(), nbrcomp + 1);
      } while (Exists(name) && type == GetType(name));
   } else
      {
	 return FAILURE;
      }
   minsize_out = min;
   number_out = nbrcomp;
   type_out = type;
   return SUCCESS;
}

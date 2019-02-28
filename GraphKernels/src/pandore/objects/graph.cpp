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
 * @author Regis Clouard - 1997-10-27
 * @author François Angot - 1999-10-11
 * @author Alexandre Duret-Lutz - 1999-10-11
 * @author Régis Clouard - 2001-04-11 (version 3.00)
 * @author Régis Clouard - 2002-12-09 (version 4.00)
 * @author Régis Clouard - 2008-01-14 (change weight type from float to double)
 * @author Régis Clouard - 2008-02-13 (add directed and undirected property)
 * @author François-Xavier Dupé - 2008-03-05 (add merge for directed graph)
 * @author François-Xavier Dupé - 2009-01-15 (extend graph representation)
 */

#include <pandore.h>
using namespace pandore;

/**
 * @file graph.cpp
 *
 */

/*
 * Adds a the node `n' in the list of adjacent nodes.
 * -> New dge with weight `w'.
 */
template< class Point >
GEdge *GNode< Point >::Add( Long n, Double w ) {
   return ( adjacents=new GEdge( n, adjacents, w ) );
}

/*
 * Adds a the node `n' in the list of adjacent nodes using specific edge.
 * -> New dge with weight `w'.
 */
template< class Point >
GEdge *GNode< Point >::Add( Long n, Long i, Double w ) {
   return ( adjacents=new GEdge( n, adjacents, i, w ) );
}

/*
 * Deletes n from the list of adjacents nodes.
 * Del does not really delete GEdge because some
 * loop variables may point to the neighbours()
 * list (ex: for (p=Neighbours(); p!=NULL; p=p->next) ...
 * @param n the node to be deleted.
 * @return	the new list of adjencent nodes.
 */
template< class Point >
GEdge *GNode< Point >::Del( Long n ) {
   GEdge *p = Neighbours();
   GEdge *q;
   
   if (!p) return NULL;

   while (p && ( p->Node() == n ))
      {
	 etrash.push_back(p);
	 p = p->Next();
	 adjacents = p;
      }

   for ( q=p; q != NULL ; q=q->Next() ) {
      if ( q->Node() == n ) {
	 p->Next(q->Next());
	 etrash.push_back(q);
	 // delete(q); // Problem when using Neighbours() list: cannot delete
      }
      else p=q;
   }
   
//    if (p->Node()==n) {
//       q=p->Next();
//       // delete(p); // Problem when using Neighbours() list: cannot delete
//       adjacents=q;
//    } else {
//       for ( q=p->Next(); q != NULL ; q=q->Next() ) {
// 	 if ( q->Node() == n ) {
// 	    p->Next(q->Next());
// 	 }
// 	 else p=q;
//       }
//    }
   return adjacents;
}

/*
 * Deletes specific edge linking to n from the list of adjacents nodes.
 * Del does not really delete GEdge because some
 * loop variables may point to the neighbours()
 * list (ex: for (p=Neighbours(); p!=NULL; p=p->next) ...
 * @param n the node to be deleted.
 * @param i the identifier of the edge
 * @return	the new list of adjencent nodes.
 */
template< class Point >
GEdge *GNode< Point >::Del( Long n, Long i ) {
   GEdge *p = Neighbours();
   GEdge *q;
   
   if (!p) return NULL;
   
   if ((p->Node()==n)&&(p->Item()==i)) {
      q=p->Next();
      etrash.push_back(p);
      // delete(p); // Problem when using Neighbours() list: cannot delete
      adjacents=q;
   } else {
      for (q=p->Next();q != NULL; p=q,q=q->Next())
	 if (q && (q->Node() ==n) && (q->Item() == i)) {
	    p->Next(q->Next());
	    etrash.push_back(q);
	    break;
	    // delete(q); // Problem when using Neighbours() list: cannot delete
	 }
   }
   return adjacents;
}

/*
 * If the node `n' is in the adjacent list
 * returns its address.
 * @param n the node number.
 */
template< class Point >
GEdge *GNode< Point >::Search( Long n ) const {
   GEdge *p;
   for (p=adjacents; p && (p->Node()!=n); p=p->Next()) ;
   return p;
}

/*
 * If the node `n' is in the adjacent list
 * returns its address.
 * @param n the node number.
 */
template< class Point >
GEdge *GNode< Point >::Search( Long n, Long i ) const {
   GEdge *p;
   for (p=adjacents; p != NULL; p=p->Next())
      if (( p->Item() == i) && (p->Node() == n))
	 break;
   return p;
}

/*
 * Node destructor.
 * -> Destroy the list of adjacent node,
 * without worrying about linked nodes.
 */
template< class Point >
GNode< Point >::~GNode( ) {
   GEdge *q,*p=adjacents;

   while ((q=p)) {
      p=p->Next();
      delete q;
   }

   std::vector<GEdge *>::iterator ite;
   for ( ite = etrash.begin(); ite != etrash.end(); ++ite )
      delete *ite;

}

/*
 * GRAPH 2D
 */

/*
 * Initializor.
 * Need a value for tnode.
 */
void Graph2d::New( Long s, Long nr, Long nc ) {
   if (tnode)
      Delete();

   ncol=nc;nrow=nr;
   if ((size=s)>0) {
      tnode=new GNode< Point2d >*[s];
      for (--s;s>=0;tnode[s--]=NULL) ;
   }
}

/*
 * Destructor.
 * Delete all nodes.
 */
void Graph2d::Delete( ) {
   if (size>0 && tnode) {
      for (int i=0;i<size;i++) {
	 if (tnode[i])
	    delete tnode[i];
      }
      delete [] tnode;
   }
   tnode=NULL;
}

/*
 * Adds a node in graph at position `node' in graph array,
 * with item=`item', and seed coordinates=`p'.
 */
Errc Graph2d::Add( Long node, Long item, const Point2d &p ) {
   if (node<size && node>=0) {
      tnode[node] = new GNode<Point2d>(item,p);
      return SUCCESS;
   }
   return FAILURE;
}

/*
 * Removes the node `n' in graph.
 * -> Remove all references in its neighbourhood nodes.
 */
Errc Graph2d::Del( Long n ) {
   GEdge *p;
   
   if (!tnode[n])
      return FAILURE;
   
   for (p=tnode[n]->Neighbours(); p!=NULL;p=p->Next()) {
      tnode[p->Node()]->Del(n);
   }
   delete tnode[n];
   tnode[n]=NULL;

   return SUCCESS;
}

/**
 * Sets the values of the graph from the
 * specified graph values.
 * @param src	the reference graph.
 */
Graph2d &Graph2d::operator=( const Graph2d &g ) {
   Long	i;
   GEdge *l;
   
   New(g.size,g.nrow,g.ncol);
   for (i=0; i<size; i++) {
      if (g[i] != NULL) {
	 Add(i,g[i]->Item(),g[i]->seed);
	 tnode[i]->value=g[i]->value;
	 for (l=g[i]->Neighbours(); l!=NULL; l=l->Next()) {
	    if (i>l->Node())
	       Link(i,l->Node(),l->Item(),l->weight);
	 }
      }
   }
   return *this;
}

/*
 * Creates a copy.
 */
Pobject *Graph2d::Clone( ) const {
   Graph2d *tmp = new Graph2d(size,nrow,ncol,_directed);
   *tmp=*this;
   return tmp;
}

/*
 * Builds a graph from a region map.
 * Do not use the region 0.
 * Use the 2nd region map to set seed cordinates.
 */
Errc Graph2d::Init( const Reg2d &rgs, const Reg2d &grm ) {
   int x,y;

   if (Init(rgs)!=SUCCESS)
      return FAILURE;

   for (y=1; y<rgs.Height()-1; y++)
      for (x=1; x<rgs.Width()-1; x++) {
	 if ((grm[y][x])) {
	    tnode[rgs[y][x]]->seed.y=y;
	    tnode[rgs[y][x]]->seed.x=x;
	 }
      }
   return SUCCESS;
}

/*
 * Builds a graph from a region map.
 * Do not use the region 0.
 * Seed on left upper coordinate.
 */
Errc Graph2d::Init( const Reg2d &rgs ) {
   register int x,y;
   Ulong a,b;

   // label + 1 to add region 0.
   New(rgs.Labels()+1,rgs.Height(),rgs.Width());

   // 1st point
   if ((a=rgs[0][0])!=0) {
      if (tnode[a]==NULL) {
	 if (!Add(a,a)) {
	    std::cerr<<"Error: Cannot allocate new node."<<std::endl;
	    return FAILURE;
	 }
	 tnode[a]->seed.x=0;
	 tnode[a]->seed.y=0;
      }
   }
   // 1st line
   for (x=1; x<rgs.Width(); x++) {
      if ((a=rgs[0][x])!=0) {
	 if (tnode[a]==NULL) {
	    if (!Add(a,a)) {
	       std::cerr<<"Error: Cannot allocate new node."<<std::endl;
	       return FAILURE;
	    }
	    tnode[a]->seed.x=x;
	    tnode[a]->seed.y=0;
	 }
	 b=rgs[0][x-1];
	 if (b && b != a) {
	    Link(a,b);
	 }
      }
   }

   for (y=1; y<rgs.Height(); y++) {
      for (x=0; x<rgs.Width(); x++) {
	 if ((a=rgs[y][x])) {
	    if (tnode[a]==NULL) {
	       if (!Add(a,a)) {
		  std::cerr<<"Error: Cannot allocate new node."<<std::endl;
		  return FAILURE;
	       }
	       tnode[a]->seed.x=x;
	       tnode[a]->seed.y=y;
	    }
	    for (int v=0; v<4; v++) { 	 // Look for 8-connected regions.
	       if ( rgs.Hold(y+v8y[v],x+v8x[v]) && ((b=rgs[y+v8y[v]][x+v8x[v]])!=0) && (b!=a)) {
		  Link(a,b);
	       }
	    }
	 }
      }
   }

   return SUCCESS;
}

/*
 * Creates a new edge between n1 and n2. If the graph
 * is undirected creates also an edge between n2 and n1.
 * If the link already exists,
 * just updates the weight value.
 * @param add if false set the value, if true add the value.
 */
Errc Graph2d::Link( Long n1, Long n2, Double w, bool add ) {
   GEdge *p;

   // New edge between n1 and n2.
   if (!(p=tnode[n1]->Search(n2)))
      tnode[n1]->Add(n2,w);
   else
      if (add) p->weight+=w; else p->weight=w;
   if ((!_directed) && (n1 != n2)){   
      // New edge between n2 and n1
      if (!(p=tnode[n2]->Search(n1)))
         tnode[n2]->Add(n1,w);
      else
	 if (add) p->weight+=w; else p->weight=w;
   }
   return SUCCESS;
}

/*
 * Creates a new edge between n1 and n2. If the graph
 * is undirected creates also an edge between n2 and n1.
 * If the link already exists,
 * just updates the weight value.
 * The edge is identified by i.
 * @param add if false set the value, if true add the value.
 */
Errc Graph2d::Link( Long n1, Long n2, Long i, Double w, bool add ) {
   GEdge *p;

   // New edge between n1 and n2.
   if (!(p=tnode[n1]->Search(n2,i)))
      tnode[n1]->Add(n2,i,w);
   else
      if (add) p->weight+=w; else p->weight=w;
   if ((!_directed) && (n1 != n2)){   
      // New edge between n2 and n1
      if (!(p=tnode[n2]->Search(n1)))
         tnode[n2]->Add(n1,i,w);
      else
	 if (add) p->weight+=w; else p->weight=w;
   }
   return SUCCESS;
}

/*
 * Removes the edges between n1 and n2,
 * and vice versa if the graph is undirected.
 */
Errc Graph2d::Unlink( Long n1, Long n2 ) {
   tnode[n1]->Del(n2);
   if ((!_directed) && (n1 != n2)) tnode[n2]->Del(n1);
   return SUCCESS;
}

/*
 * Removes the specified edge between n1 and n2,
 * and vice versa if the graph is undirected.
 */
Errc Graph2d::Unlink( Long n1, Long n2, Long i ) {
   tnode[n1]->Del(n2,i);
   if ((!_directed) && (n1 != n2)) tnode[n2]->Del(n1,i);
   return SUCCESS;
}

/*
 * Merges 2 nodes in 1 ->n1
 * Updates the adjacent list of all node (n2 ->n1).
 * and delete n2.
 */
Errc Graph2d::Merge( Long n1, Long n2 ) {
   GEdge *ptr;
   
   if (_directed) {
      // std::cerr << "Error: Merge is not implemented for directed graphs." << std::endl;
      // return FAILURE;
      
      if ((n1==n2) || (!tnode[n1]) || (!tnode[n2])) return FAILURE;

      // Remove links
      Unlink(n1,n2);
      Unlink(n2,n1);

      while((ptr=tnode[n2]->Neighbours())){
	 if ( n2 != ptr->Node() )
	    Link(n1,ptr->Node(),ptr->Item(),ptr->weight,true);
	 else
	    Link(n1,n1,ptr->Item(),ptr->weight,true);
	 Unlink(n2,ptr->Node(),ptr->Item());
      }

      // Cross the graph in order to replace n2 by n1
      for ( int i = 0; i < this->Size(); ++i ) {
	 const GNode<Point2d> * noeud = (*this)[i];
	 if ( noeud != NULL ) {
	    for ( ptr = noeud->Neighbours() ; ptr != NULL; ptr = ptr->Next() )
	       if ( ptr->Node() == n2 ) {
		  Link(i,n1,ptr->Item(),ptr->weight,true);
		  Unlink(i,n2,ptr->Item());
	       }
	 }
      }
      Del(n2);
   } else {
      if ((n1==n2) || (!tnode[n1]) || (!tnode[n2])) return FAILURE;
      
      Unlink(n1,n2);
      
      while((ptr=tnode[n2]->Neighbours())){
	 if ( n2 != ptr->Node() )
	    Link(n1,ptr->Node(),ptr->Item(),ptr->weight,true);
	 else
	    Link(n1,n1,ptr->Item(),ptr->weight,true);
	 Unlink(n2,ptr->Node(),ptr->Item());
      }
      Del(n2);
   }
   return SUCCESS;
}

/*
 * Splits the node n1 in n2.
 * Updates all adjacent lists.
 */
Errc Graph2d::Split( Long n1, Long n2 ) {
   Add(n2,tnode[n1]->Item(),tnode[n1]->seed);
   for (GEdge *ptr=tnode[n1]->Neighbours();ptr;ptr->Next())
      if ( n1 != ptr->Node() )
	 Link(n2,ptr->Node(),ptr->Item());
      else
	 Link(n2,n2,ptr->Item(),ptr->weight);
   return SUCCESS;
}

/*      
 * Loads graph atributes from the file.
 */
Errc Graph2d::LoadAttributes( FILE *df ) {
   Long attr[3];
   char d;

   if (Fdecode(attr,sizeof(Long),3,df) < 3)
      return FAILURE;
   if (Fdecode(&d,sizeof(d),1,df) < 1)
      return FAILURE;
   New(attr[0],attr[1],attr[2]);
   _directed = (d!=0);
   return SUCCESS;
}

/*
 * Saves graph attributes into the file.
 */
Errc Graph2d::SaveAttributes( FILE *df ) const {
   Long attr[3];
   char d= _directed;

   attr[0]=size; attr[1]=nrow;attr[2]=ncol;
   if (Fencode(attr,sizeof(Long),3,df) < 3)
      return FAILURE;
   if (Fencode(&d,sizeof(d),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Loads a graph from a file.
 */
Errc Graph2d::LoadData( FILE *df ) {
   register int i;
   Long suiv;
   Long item;
   Graph2d::ValueType val;
   Graph2d::ValueType wgh;
   Point2d p;

   for (i=0;i<size;i++) {
      if (Fdecode(&suiv,sizeof(suiv),1,df) < 1)
	 return FAILURE;
      if (suiv >= 0) {
	 if (Fdecode(&val,sizeof(val),1,df)<1)
	    return FAILURE;
	 if (p.Load(df,_inversionMode) == FAILURE)
	    return FAILURE;
	 if (!Add(i,suiv,p)) return FAILURE;
	 tnode[i]->value=val;
	 if (Fdecode(&suiv,sizeof(suiv),1,df)<1) {
	    return FAILURE;
	 }
	 while (suiv>-1){
	    if (Fdecode(&wgh,sizeof(wgh),1,df)<1) {
	       return FAILURE;
	    }
	    if (Fdecode(&item,sizeof(item),1,df)<1)
	       return FAILURE;
	    if (i>=suiv || _directed)	// Link only with former nodes for directed graphs.
	       Link(i,suiv,item,wgh);
	    if (Fdecode(&suiv,sizeof(suiv),1,df)<1)
	       return FAILURE;
	 }
      } else tnode[i]=NULL;
   }
   
   return SUCCESS;
}

/*
 * Saves graph data into the file.
 */
Errc Graph2d::SaveData( FILE *df ) const {
   register int i;
   register GEdge *v;
   //   const Long fin=-1;
   const Long fin=-1;
   Long suiv;
   Long item;

   for (i=0;i<size;i++) {
      if (tnode[i]==NULL) {
	 if (Fencode(&fin,sizeof(fin),1,df)<1) {
	    return FAILURE;
	 }
      } else {
	 suiv=tnode[i]->Item();
	 if (Fencode(&suiv,sizeof(suiv),1,df)<1)
	    return FAILURE;
	 if (Fencode(&(tnode[i]->value),sizeof(tnode[i]->value),1,df)<1)
	    return FAILURE;
	 if (tnode[i]->seed.Save(df) == FAILURE)
	    return FAILURE;
	 for (v=tnode[i]->Neighbours();v!=NULL;v=v->Next()){
	    suiv=v->Node();
	    item=v->Item();
	    if (Fencode(&suiv,sizeof(suiv),1,df)<1)
	       return FAILURE;
	    if (Fencode(&(v->weight),sizeof(v->weight),1,df)<1)
	       return FAILURE;
	    if (Fencode(&item,sizeof(item),1,df)<1)
	       return FAILURE;
	 }
	 if (Fencode(&fin,sizeof(fin),1,df)<1)
	    return FAILURE;
      }
   }

   return SUCCESS;
}

/*
 * Creates a new graph only with the unmasked nodes by `mask'.
 */
Pobject *Graph2d::Mask( const Pobject *mask ) {
   if ((!mask)||(mask->Type()!=Po_Reg2d)||(((Reg2d*)mask)->Size()!=ImageSize())) {
      std::cerr << "Warning: bad mask format... ignored" << std::endl;
      return this;
   }

   Graph2d *objd = (Graph2d*)Clone();
   Reg2d *m=(Reg2d*)mask;
   bool *regions=(bool*)calloc(m->Labels()+1,sizeof(bool));
   Ulong *pm=m->Vector();
   
   for (int i=0; i<nrow*ncol; i++, pm++)
      regions[*pm]=true; // List of masked/unmasked regions.

   for (int g=1;g<objd->Size();g++)
      if (!regions[g]) // Delete masked nodes.
	 objd->Del(g);

   delete [] regions;
   return objd;
}

/*
 * Creates the output graph with the previous masked node 
 * and the new nodes.
 */
Pobject *Graph2d::UnMask( const Pobject *mask, const Pobject *gri ) {
   if ((!mask)||(mask->Type()!=Po_Reg2d)||(((Reg2d*)mask)->Size()!=ImageSize())||(gri->Type() != Type())||(((Graph2d*)gri)->ImageSize()!=ImageSize())){
      std::cerr << "Warning: bad unmask format... ignored" << std::endl;
      return this;
   }
   if ((mask == NULL) || (mask->Type() != Po_Reg2d) || (gri->Type() != Type())) {
      return this;
   }

   Graph2d *objs = (Graph2d*)gri;
   Graph2d *grs = (Graph2d*)this;
   Graph2d *objd = (Graph2d*)gri->Clone();
   Reg2d *m=(Reg2d*)mask;
   bool *regions=(bool*)calloc(m->Labels()+1,sizeof(bool));
   Ulong *pm=m->Vector();
   
   for (int i=0; i<nrow*ncol; i++, pm++)
      regions[*pm]=true; // List of masked/unmasked regions.

   // Add or remove node from the input graph with processed graph.
   for (int g=1;g<objd->Size();g++){
      if ( !((*objs)[g]) && (*grs)[g]){  // Added node
	 objd->Add(g,(*grs)[g]->Item(),(*grs)[g]->seed);
	 for (GEdge *ptr=(*grs)[g]-> Neighbours(); ptr!=NULL; ptr=ptr->Next())
	    objd->Link(g,ptr->Node(),ptr->Item(),ptr->weight);
      }
      if (regions[g]){
	 if ( (*objs)[g] && !((*grs)[g]) ) // Removed node.
	    objd->Del(g);
	 if ((*objs)[g] && (*grs)[g] ) // Update neighbours.
	    for (GEdge *ptr=(*grs)[g]-> Neighbours(); ptr!=NULL; ptr=ptr->Next())
	       objd->Link(g,ptr->Node(),ptr->Item(),ptr->weight);
      }
   }
   delete [] regions;
   return objd;
}


/*
 * GRAPH 3D 
 */

/*
 * Initializator.
 * Need a value for tnode.
 */
void Graph3d::New( Long s, Long nd, Long nr, Long nc ) {
   if (tnode)
      Delete();
   
   ncol=nc;nrow=nr;ndep=nd;
   if ((size=s)>0){
      tnode=new GNode< Point3d >*[s];
      for (--s;s>=0;tnode[s--]=NULL) ;
   }
}

/*
 * Destructor.
 * Delete all nodes.
 */
void Graph3d::Delete( ) {
   if (size>0 && tnode){
      for (int i=0;i<size;i++){
	 if (tnode[i])
	    delete tnode[i];
      }
      delete [] tnode;
   }
   tnode=NULL;
}

/*
 * Enlarge the number of possible nodes in the graph
 */
Errc Graph3d::Enlarge( Long s){
   if(s > size){
      GNode<Point3d> ** tnode_tmp = new GNode<Point3d> * [s];
      memset(tnode_tmp,0,sizeof(GNode<Point3d>*)*s);
      for(int i=0;i<size;i++)
	 tnode_tmp[i] = tnode[i];
      delete [] tnode;
      tnode = tnode_tmp;
      size = s;
   }
}
/*
 * Adds a node in graph
 * at position `node' in graph array,
 * with item=`item',
 * and seed coordinates=`p'.
 */
Errc Graph3d::Add( Long node, Long item, const Point3d &p ) {
   if (node<size && node>=0){
      tnode[node] = new GNode<Point3d>(item,p);
      return SUCCESS;
   }
   return FAILURE;
}

/*
 * Removes the node `n' in graph.
 * -> Removes all refences in its
 * neighbourhood nodes.
 */
Errc Graph3d::Del( Long n ) {
   GEdge *p;
   
   if (!tnode[n])
      return FAILURE;
   
   for (p=tnode[n]->Neighbours(); p!=NULL;p=p->Next()){
      tnode[p->Node()]->Del(n);
   }
   delete tnode[n];
   tnode[n]=NULL;

   return SUCCESS;
}

/*
 * Creates a copy of graph g.
 */
Graph3d &Graph3d::operator=( const Graph3d &g ) {
   Long	i;
   GEdge *l;

   New(g.size,g.ndep,g.nrow,g.ncol);
   for (i=0; i<size; i++){
      if (g[i] != NULL){
	 Add(i,g[i]->Item(),g[i]->seed);
	 tnode[i]->value=g[i]->value;
	 for (l=g[i]->Neighbours(); l!=NULL; l=l->Next()){
	    if (i>l->Node())
	       Link(i,l->Node(),l->Item(),l->weight);
	 }
      }
   }
   return *this;
}

/*
 * Creates a copy.
 */
Pobject *Graph3d::Clone( ) const {
   Graph3d *tmp = new Graph3d(size,ndep,nrow,ncol,_directed);
   *tmp=*this;
   return tmp;
}

/*
 * Builds a graph from a region map.
 * Do not use the region 0.
 * Use the 2nd region map as seed cordinates.
 */
Errc Graph3d::Init( const Reg3d &rgs, const Reg3d &grm ) {
   int x,y,z;
   
   if (Init(rgs) != SUCCESS)
      return FAILURE;
   
   for (z=1; z<rgs.Depth()-1; z++)
      for (y=1; y<rgs.Height()-1; y++)
	 for (x=1; x<rgs.Width()-1; x++){
	    if ((grm[z][y][x])){
	       tnode[rgs[z][y][x]]->seed.z=z;
	       tnode[rgs[z][y][x]]->seed.y=y;
	       tnode[rgs[z][y][x]]->seed.x=x;
	    }
	 }
   
   return SUCCESS;
}

/*
 * Builds a graph from a region map.
 * Do not use the region 0.
 * No seed.
 */
Errc Graph3d::Init( const Reg3d &rgs ) {
   register short int x,y,z;
   Ulong a,b;
   
   // labels + 1 to add the region 0
   New(rgs.Labels()+1,rgs.Depth(),rgs.Height(),rgs.Width());
   
   // 1st point
   if ((a=rgs[0][0][0])!=0) {
      if (tnode[a]==NULL) {
	 if (!Add(a,a)) {
	    std::cerr<<"Error: Cannot allocate new node."<<std::endl;
	    return FAILURE;
	 }
	 tnode[a]->seed.x=0;
	 tnode[a]->seed.y=0;
	 tnode[a]->seed.z=0;
      }
   }

   // 1st line
   for (x=1; x<rgs.Width(); x++){
      if ((a=rgs[0][0][x])!=0) {
	 if (tnode[a]==NULL) {
	    if (!Add(a,a)) {
	       std::cerr<<"Error: Cannot allocate new node."<<std::endl;
	       return FAILURE;
	    }
	    tnode[a]->seed.x=x;
	    tnode[a]->seed.y=0;
	    tnode[a]->seed.z=0;
	 }
	 if ((b=rgs[0][0][x-1]) != a) {
	    Link(a,b);
	 }
      }
   }
   // 1st plane.
   for (y=1; y<rgs.Height()-1; y++){
      for (x=0; x<rgs.Width()-1; x++) {
	 if ((a=rgs[0][y][x])) {
	    if (tnode[a]==NULL) {
	       if (!Add(a,a)) {
		  std::cerr<<"Error: Cannot allocate new node."<<std::endl;
		  return FAILURE;
	       }
	       tnode[a]->seed.x=x;
	       tnode[a]->seed.y=y;
	       tnode[a]->seed.z=0;
	    }
	    for (int v=9; v<13; v++) { 	 // Look for 8-connected regions.
	       if ( rgs.Hold(0,y+v26y[v],x+v26x[v]) && ((b=rgs[0][y+v26y[v]][x+v26x[v]])!=0) && (b!=a)) {
		  Link(a,b);
	       }
	    }
	 }
      }
   }
   // general case.
   for (z=1; z<rgs.Depth(); z++) {
      for (y=0; y<rgs.Height(); y++) {
	 for (x=0; x<rgs.Width(); x++) {
	    if ((a=rgs[z][y][x])) {
	       if (tnode[a]==NULL) {
		  if (!Add(a,a)) {
		     std::cerr<<"Error: Cannot allocate new node."<<std::endl;
		     return FAILURE;
		  }
		  tnode[a]->seed.x=x;
		  tnode[a]->seed.y=y;
		  tnode[a]->seed.z=z;
	       }
	       for (int v=0; v<13; v++) { 	 // Look for 8-connected regions.
		  if ( rgs.Hold(z+v26z[v],y+v26y[v],x+v26x[v])
		       && ((b=rgs[z+v26z[v]][y+v26y[v]][x+v26x[v]])!=0) && (b!=a)) {
		     Link(a,b);
		  }
	       }
	    }
	 }
      }
   }
   return SUCCESS;
}

/*
 * Creates a new edge between n1 and n2. If the graph
 * is undirected creates also an edge between n2 and n1.
 * If the link already exists,
 * just updates the weight value.
 */
Errc Graph3d::Link( Long n1, Long n2, Double w, bool add ) {
   GEdge *p;

   if (!(p=tnode[n1]->Search(n2)))
      tnode[n1]->Add(n2,w);
   else
      if (add) p->weight+=w; else p->weight=w;

   if ((!_directed ) && (n1 != n2)) {   
      if (!(p=tnode[n2]->Search(n1)))
	 tnode[n2]->Add(n1,w);
      else
	 if (add) p->weight+=w; else p->weight=w;
   }
   return SUCCESS;
}

/*
 * Creates a new edge between n1 and n2. If the graph
 * is undirected creates also an edge between n2 and n1.
 * If the link already exists,
 * just updates the weight value.
 * The edge has a specific identifier.
 */
Errc Graph3d::Link( Long n1, Long n2, Long i, Double w, bool add ) {
   GEdge *p;
   
   if (!(p=tnode[n1]->Search(n2,i)))
      tnode[n1]->Add(n2,i,w);
   else
      if (add) p->weight+=w; else p->weight=w;

   if ((!_directed ) && (n1 != n2)) {   
      if (!(p=tnode[n2]->Search(n1,i)))
	 tnode[n2]->Add(n1,i,w);
      else
	 if (add) p->weight+=w; else p->weight=w;
   }
   return SUCCESS;
}

/*
 * Removes the edge between n1 and n2
 * and vice versa if the graph is undirected.
 */
Errc Graph3d::Unlink( Long n1, Long n2 ) {
   tnode[n1]->Del(n2);
   if ((!_directed) && (n1 != n2)) tnode[n2]->Del(n1);
   return SUCCESS;
}

/*
 * Removes the identified edge between n1 and n2
 * and vice versa if the graph is undirected.
 */
Errc Graph3d::Unlink( Long n1, Long n2, Long i ) {
   tnode[n1]->Del(n2,i);
   if ((!_directed) && (n1 != n2)) tnode[n2]->Del(n1,i);
   return SUCCESS;
}

/*
 * Merges 2 nodes in 1 ->n1
 * Updates the adjacent list of all node (n2 ->n1).
 * and delete n2.
 */
Errc Graph3d::Merge( Long n1, Long n2 ) {
   GEdge *ptr;

   if (_directed) {
       if ((n1==n2) || (!tnode[n1]) || (!tnode[n2])) return FAILURE;

      // Remove links
      Unlink(n1,n2);
      Unlink(n2,n1);

      while((ptr=tnode[n2]->Neighbours())){
	 if ( n2 != ptr->Node() )
	    Link(n1,ptr->Node(),ptr->Item(),ptr->weight,true);
	 else
	    Link(n1,n1,ptr->Item(),ptr->weight,true);
	 Unlink(n2,ptr->Node(),ptr->Item());
      }

      // Cross the graph in order to replace n2 by n1
      for ( int i = 0; i < this->Size(); ++i ) {
	 const GNode<Point3d> * noeud = (*this)[i];
	 if ( noeud != NULL ) {
	    for ( ptr = noeud->Neighbours() ; ptr != NULL; ptr = ptr->Next() )
	       if ( ptr->Node() == n2 ) {
		  Link(i,n1,ptr->Item(),ptr->weight,true);
		  Unlink(i,n2,ptr->Item());
	       }
	 }
      }
      Del(n2);
   }  else {
      if ((n1==n2) || (!tnode[n1]) || (!tnode[n2])) return FAILURE;
      
      Unlink(n1,n2);
      
      while((ptr=tnode[n2]->Neighbours())){
	 if ( n2 != ptr->Node() )
	    Link(n1,ptr->Node(),ptr->Item(),ptr->weight,true);
	 else
	    Link(n1,n1,ptr->Item(),ptr->weight,true);
	 Unlink(n2,ptr->Node(),ptr->Item());
      }
      Del(n2);
   }

   return SUCCESS;
}

/*
 * Splits the node n1 in n2.
 * Updates all adjacent lists.
 */
Errc Graph3d::Split( Long n1, Long n2 ) {
   Add(n2,tnode[n1]->Item(),tnode[n1]->seed);
   for (GEdge *ptr=tnode[n1]->Neighbours();ptr;ptr->Next())
      if ( n1 != ptr->Node() )
	 Link(n2,ptr->Node(),ptr->Item());
      else
	 Link(n2,n2,ptr->Item(),ptr->weight);
   return SUCCESS;
}

/*      
 * Loads graph atributes from the file.
 */
Errc Graph3d::LoadAttributes( FILE *df ) {
   Long attr[4];
   bool d;
   
   if (Fdecode(attr,sizeof(Long),4,df) < 4)
      return FAILURE;
   if (Fdecode(&d,sizeof(bool),1,df) < 1)
      return FAILURE;
   _directed = d;
   New(attr[0],attr[1],attr[2],attr[3]);
   
   return SUCCESS;
}

/*
 * Saves graph attributes into the file.
 */
Errc Graph3d::SaveAttributes( FILE *df ) const {
   Long attr[4];
   bool d = _directed;
   attr[0]=size; attr[1]=ndep; attr[2]=nrow;attr[3]=ncol;
   if (Fencode(attr,sizeof(Long),4,df) < 4)
      return FAILURE;
   if (Fencode(&d,sizeof(bool),1,df) < 1)
      return FAILURE;
   return SUCCESS;
}

/*
 * Loads a graph from a file.
 */
Errc Graph3d::LoadData( FILE *df ) {
   register int i;
   Long suiv;
   Long item;
   Graph2d::ValueType val;
   Graph2d::ValueType wgh;
   Point3d p;

   for (i=0;i<size;i++){
      if (Fdecode(&suiv,sizeof(suiv),1,df) < 1)
	 return FAILURE;
      if (suiv >= 0){
	 if (Fdecode(&val,sizeof(val),1,df)<1)
	    return FAILURE;
	 if (p.Load(df,_inversionMode) == FAILURE)
	    return FAILURE;
	 if (!Add(i,suiv,p)) return FAILURE;
	 tnode[i]->value=val;
	 
	 if (Fdecode(&suiv,sizeof(suiv),1,df)<1)
	    return FAILURE;
	 while (suiv>-1){
	    if (Fdecode(&wgh,sizeof(wgh),1,df)<1)
	       return FAILURE;
	    if (Fdecode(&item,sizeof(item),1,df)<1)
	       return FAILURE;
	    if (i>=suiv || _directed)	// Link only with former nodes for directed graphs.
	       Link(i,suiv,item,wgh);
	    if (Fdecode(&suiv,sizeof(suiv),1,df)<1)
	       return FAILURE;
	 }
      } else tnode[i]=NULL;
   }
   
   return SUCCESS;
}

/*
 * Saves graph data into the file.
 */
Errc Graph3d::SaveData( FILE *df ) const {
   register int i;
   register GEdge *v;
   Long fin=-1;
   Long suiv;
   Long item;

   for (i=0;i<size;i++) {
      if (tnode[i]==NULL) {
	 if (Fencode(&fin,sizeof(fin),1,df)<1)
	    return FAILURE;
      } else {
	 suiv=tnode[i]->Item();
	 if (Fencode(&suiv,sizeof(suiv),1,df)<1)
	    return FAILURE;
	 if (Fencode(&(tnode[i]->value),sizeof(tnode[i]->value),1,df)<1)
	    return FAILURE;
	 if (tnode[i]->seed.Save(df) == FAILURE)
	    return FAILURE;
	 for (v=tnode[i]->Neighbours();v!=NULL;v=v->Next()) {
	    suiv=v->Node();
	    item=v->Item();
	    if (Fencode(&suiv,sizeof(suiv),1,df)<1)
	       return FAILURE;
	    if (Fencode(&(v->weight),sizeof(v->weight),1,df)<1)
	       return FAILURE;
	    if (Fencode(&item,sizeof(item),1,df)<1)
	       return FAILURE;
	 }
	 if (Fencode(&fin,sizeof(fin),1,df)<1)
	    return FAILURE ;
      }
   }
   return SUCCESS;
}

/*
 * Creates a new graph only with the unmasked nodes by `mask'.
 */
Pobject *Graph3d::Mask( const Pobject *mask ) {
   if ((!mask)||(mask->Type()!=Po_Reg3d)||(((Reg3d*)mask)->Size()!=ImageSize())){
      std::cerr << "Warning: bad mask format... ignored" << std::endl;
      return this;
   }

   Graph3d *objd = (Graph3d*)Clone();
   Reg3d *m=(Reg3d*)mask;
   bool *regions=(bool*)calloc(m->Labels()+1,sizeof(bool));
   Ulong *pm=m->Vector();
   
   for (int i=0; i<ndep*nrow*ncol; i++, pm++)
      regions[*pm]=true; // List of masked/unmasked regions.

   for (int g=1;g<objd->Size();g++)
      if (!regions[g]) // Delete masked nodes.
	 objd->Del(g);

   delete [] regions;
   return objd;
}

/*
 * CreateS the output graph with the previous masked node 
 * and the new nodes.
 */
Pobject *Graph3d::UnMask( const Pobject *mask, const Pobject *gr ) {
   if ((!mask)||(mask->Type()!=Po_Reg3d)||(((Reg3d*)mask)->Size()!=ImageSize())||(gr->Type() != Type())||(((Graph3d*)gr)->ImageSize()!=ImageSize())){
      std::cerr << "Warning: bad unmask format... ignored" << std::endl;
      return this;
   }
   if ((mask == NULL) || (mask->Type() != Po_Reg3d) || (gr->Type() != Type())){
      return this;
   }

   Graph3d *objs = (Graph3d*)gr;
   Graph3d *grs = (Graph3d*)this;
   Graph3d *objd = (Graph3d*)gr->Clone();
   Reg3d *m=(Reg3d*)mask;
   bool *regions=(bool*)calloc(m->Labels()+1,sizeof(bool));
   Ulong *pm=m->Vector();
   
   for (int i=0; i<ndep*nrow*ncol; i++, pm++)
      regions[*pm]=true; // List of masked/unmasked regions.

   // Add or remove node from the input graph with processed graph.
   for (int g=1;g<objd->Size();g++){
      if ( !((*objs)[g]) && (*grs)[g]){  // Added node
	 objd->Add(g,(*grs)[g]->Item(),(*grs)[g]->seed);
	 for (GEdge *ptr=(*grs)[g]-> Neighbours(); ptr!=NULL; ptr=ptr->Next())
	    objd->Link(g,ptr->Node(),ptr->Item(),ptr->weight);
      }
      if (regions[g]){
	 if ( (*objs)[g] && !((*grs)[g]) ) // Removed node.
	    objd->Del(g);
	 if ((*objs)[g] && (*grs)[g] ) // Update neighbours.
	    for (GEdge *ptr=(*grs)[g]-> Neighbours(); ptr!=NULL; ptr=ptr->Next())
	       objd->Link(g,ptr->Node(),ptr->Item(),ptr->weight);
      }
   }
   delete [] regions;
   return objd;
}


// Explicit instanciation: not necessary with Visual
template class pandore::GNode< Point2d >;
template class pandore::GNode< Point3d >;

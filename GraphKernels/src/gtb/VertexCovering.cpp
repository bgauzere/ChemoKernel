/*-----------------------------------------------------------------------

  File        : VertexCovering.cpp

  Description : Vertex covering by paths on trees

  This is based on the paper "Two fixed-parameter algorithms for
  Vertex Covering by Paths on Trees" by
  Jiong Guo, Rolf Niedermeier and Johannes Uhlmann, in Information
  Processing Letters 106(2008) 81-86

  Copyright  : Francois-Xavier Dupé - http://www.greyc.ensicaen.fr/~fdupe/
               	       
  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.

  ---------------------------------------------------------------------*/

#include "VertexCovering.h"
#include <queue>
#include <map>
#include <limits>

using namespace pandore;
using namespace std;
using namespace cimg_library;

void VertexCovering::getVertexSuffixTree ( Graph3d * grp, deque<int> & tab, int node, CImg<bool> & visit )
{
  visit(node) = true;
  
  GEdge * edge = (*grp)[node]->Neighbours();
  for ( ; edge != 0; edge = edge->Next() )
    {
      if ( ! visit(edge->Node()) )
	getVertexSuffixTree(grp,tab,edge->Node(),visit);
    }
  
  tab.push_back(node);
}

deque<int> VertexCovering::getVertexSuffixTree ( void )
{
  Graph3d * grp = col->GETPOBJECT("tree",Graph3d);
  CImg<bool> visit(grp->Size(),1,1,1,false);
  deque<int> res;
  int node = -1;
  int othernode = 0;
  
  // 1 - Seeking for the first non-null node
  for ( int i = 0; i < grp->Size(); ++i )
    if ( (*grp)[i] != 0 )
      {
	othernode = i; // In case of no good node

	int count = 0;
	GEdge *edge = (*grp)[i]->Neighbours();
	for ( ; edge != 0; edge = edge->Next() )
	  ++count;

	if ( count > 2 )
	  {
	    node = i;
	    break;
	  }
      }

  if ( node == -1 )
    node = othernode;

  visit(node) = true;

  // 2 - Looking at the children
  GEdge * edge = (*grp)[node]->Neighbours();
  for ( ; edge != 0; edge = edge->Next() )
    {
      if ( ! visit(edge->Node()) )
	getVertexSuffixTree(grp,res,edge->Node(),visit);
    }
  
  res.push_back(node);

  return res;
}

void VertexCovering::computeVpath ( Graph3d * grp, const deque<trail> & bag, const deque<int> & suffix, map<int,NodeInfo> & sNode )
{
  CImg<bool> visit ( grp->Size(), 1, 1, 1, false );
  
  for ( unsigned int i = 0; i < suffix.size(); ++i )
    {
      deque<trail> A; // paths covering at least one edge to a child of suffix[i]
      deque<trail> B; // paths covering (exactly) two children of suffix[i]

      visit(suffix[i]) = true;
      
      // Compute R => the set of paths with one endpoint being suffix[i]
      deque<trail> R;
      for ( unsigned int j = 0; j < bag.size(); ++j )
	{
	  if (( bag[j].t.front() == suffix[i] )||( bag[j].t.back() == suffix[i] ))
	    R.push_back(bag[j]);
	}

      // Looking at children
      GEdge * edge = (*grp)[suffix[i]]->Neighbours();
      for ( ; edge != 0; edge = edge->Next() )
	if ( visit(edge->Node()) )
	  {
	    B = unionset(B,intersection(A,sNode[edge->Node()].epath));
	    A = unionset(A,sNode[edge->Node()].epath);
	  }

      sNode[suffix[i]].vpath = unionset(A,R);

      if ( suffix[i] != suffix.back() )
	sNode[suffix[i]].epath = unionset(difference(R,A),difference(A,unionset(B,R)));
    }
}

deque<trail> VertexCovering::unionset ( const deque<trail> & s1, const deque<trail> & s2 )
{
  deque<trail> res = s1;
  
  for ( unsigned int i = 0; i < s2.size(); ++i )
    {
      bool addpath = true;

      for ( unsigned int j = 0; j < s1.size(); ++j )
	if ((( s2[i].t.front() == s1[j].t.front() ) && ( s2[i].t.back() == s1[j].t.back() )) ||
	    (( s2[i].t.back() == s1[j].t.front() ) && ( s2[i].t.front() == s1[j].t.back() )))
	  {
	    addpath = false;
	    break;
	  }

      if ( addpath )
	res.push_back(s2[i]);
    }

  return res;
}

deque<trail> VertexCovering::intersection ( const deque<trail> & s1, const deque<trail> & s2 )
{
  deque<trail> res;
  
  for ( unsigned int i = 0; i < s1.size(); ++i )
    {
      bool addpath = false;

      for ( unsigned int j = 0; j < s2.size(); ++j )
	if ((( s1[i].t.front() == s2[j].t.front() ) && ( s1[i].t.back() == s2[j].t.back() )) ||
	    (( s1[i].t.back() == s2[j].t.front() ) && ( s1[i].t.front() == s2[j].t.back() )))
	  {
	    addpath = true;
	    break;
	  }

      if ( addpath )
	res.push_back(s1[i]);
    }

  return res;
}

deque<trail> VertexCovering::difference ( const deque<trail> & s1, const deque<trail> & s2 )
{
  deque<trail> res;
  
  for ( unsigned int i = 0; i < s1.size(); ++i )
    {
      bool addpath = true;

      for ( unsigned int j = 0; j < s2.size(); ++j )
	if ((( s1[i].t.front() == s2[j].t.front() ) && ( s1[i].t.back() == s2[j].t.back() )) ||
	    (( s1[i].t.back() == s2[j].t.front() ) && ( s1[i].t.front() == s2[j].t.back() )))
	  {
	    addpath = false;
	    break;
	  }

      if ( addpath )
	res.push_back(s1[i]);
    }
  
  return res;
}

deque<trail> VertexCovering::pathVertexCovering ( const deque<int> & suffix, map<int,NodeInfo> & sNode )
{
  Graph3d * grp = col->GETPOBJECT("tree",Graph3d);
  CImg<bool> visit (grp->Size(),1,1,1,false);
 
  // Complete the whole dynamic array
  for ( unsigned int i = 0; i < suffix.size(); ++i )
    {
      visit(suffix[i]) = true;

      GEdge * edge = (*grp)[suffix[i]]->Neighbours();
      for ( ; edge != 0 ; edge = edge->Next() )
	if ( visit(edge->Node()) )
	  {
	    // Compute masks c1 and c2
	    unsigned int c1 = 0; // Path belonging to the node and the current child
	    unsigned int c2 = 0; // Path only in the current node
	    unsigned int maskbis = 0; // The common paths between child and parent for the child

	    for ( unsigned int x = 0; x < sNode[suffix[i]].vpath.size(); ++x )
	      {
		bool addinc1 = false;
		trail & s1 = sNode[suffix[i]].vpath[x];

		for ( unsigned int y = 0; y < sNode[edge->Node()].vpath.size(); ++y )
		  {
		    trail & s2 = sNode[edge->Node()].vpath[y];
		    
		    if ((( s1.t.front() == s2.t.front() ) && ( s1.t.back() == s2.t.back() )) ||
			(( s1.t.back() == s2.t.front() ) && ( s1.t.front() == s2.t.back() )))
		      {
			addinc1 = true;
			maskbis = maskbis | ( 1 << y );
			break;
		      }
		  }

		unsigned index = 1 << x;		
		if ( addinc1 )
		  c1 = c1 | index; // Add the path inside c1
		else
		  c2 = c2 | index; // Add the path inside c1
	      }

	    // Update the dynamic table
	    unsigned int size = 1 << sNode[suffix[i]].vpath.size();
	    for ( unsigned int x = 0; x < size; x++ )
	      {
		// Weight of the path in c1
		double c1sum = 0.0;
		unsigned int c1bis = 0;
		for ( unsigned int y = 0; y < sNode[suffix[i]].vpath.size(); y++ )
		  {
		    unsigned int index = 1 << y;
		    index = index & c1;
		    
		    if (( x & index ) != 0 )
		      {
			c1sum += sNode[suffix[i]].wpath[y];

			// Seeking for the equivalent path in child
			for ( unsigned int z = 0; z < sNode[edge->Node()].vpath.size(); z++ )
			  {
			    if ((( sNode[edge->Node()].vpath[z].t.front() == sNode[suffix[i]].vpath[y].t.front()) &&
				 ( sNode[edge->Node()].vpath[z].t.back()  == sNode[suffix[i]].vpath[y].t.back() )) ||
				(( sNode[edge->Node()].vpath[z].t.front() == sNode[suffix[i]].vpath[y].t.back() ) &&
				 ( sNode[edge->Node()].vpath[z].t.back()  == sNode[suffix[i]].vpath[y].t.front() )))
			      c1bis = c1bis | ( 1 << z );
			  }
		      }
		  }
		
		// Minimum weight inside combinaisons from child
		unsigned int csize = 1 << sNode[edge->Node()].vpath.size();
		double minimum = numeric_limits<double>::max();
		for ( unsigned int z = 0; z < csize; ++z )
		  {
		    unsigned int index = ( z & (~maskbis) ) | c1bis;
		    if ( sNode[edge->Node()].dynaTab[index] < minimum )
		      minimum = sNode[edge->Node()].dynaTab[index];
		  }

		// Update
		sNode[suffix[i]].dynaTab[x] += minimum - c1sum;
	      }
	  }
      
    }

  // Get the optimal covering
  unsigned int root = suffix.back();
  deque<trail> res;
  visit.fill(false);
  computeVertexSolution(root,sNode,res,visit);
  
  return res;
}

void VertexCovering::computeVertexSolution ( int node, map<int,NodeInfo> & sNode, deque<trail> & res, CImg<bool> & visit )
{
  Graph3d * grp = col->GETPOBJECT("tree",Graph3d);
  
  // Compute the current mask
  unsigned int mask = 0;
  for ( unsigned int i = 0; i < sNode[node].vpath.size(); ++i )
    {
      bool addpath = false;
      
      for ( unsigned int j = 0; j < res.size(); ++j )
	if ((( sNode[node].vpath[i].t.front() == res[j].t.front() ) &&
	     ( sNode[node].vpath[i].t.back()  == res[j].t.back() )) ||
	    (( sNode[node].vpath[i].t.back()  == res[j].t.front() ) &&
	     ( sNode[node].vpath[i].t.front() == res[j].t.back() )))
	  {
	    addpath = true;
	    break;
	  }

      if ( addpath )
	mask = mask | ( 1 << i );
    }
    
  // Find the path set which gives the minimal weight
  unsigned int size = 1 << sNode[node].vpath.size();
  unsigned int c = 0;
  double minimum = numeric_limits<double>::max();
  for ( unsigned int x = 0; x < size; ++x )
    {
      unsigned int index = x | mask;
      if ( minimum > sNode[node].dynaTab[index] )
	{
	  minimum = sNode[node].dynaTab[index];
	  c = index;
	}
    }

  // Get the corresponding path
  deque<trail> npath;
  for ( unsigned int x = 0; x < sNode[node].vpath.size(); ++x )
    if ( (( 1 << x ) & c) != 0 )
      npath.push_back(sNode[node].vpath[x]);

  // Update
  res = unionset(res,npath);

  // Look after the childs
  visit(node) = true;
  GEdge * edge = (*grp)[node]->Neighbours();
  for ( ; edge != 0; edge = edge->Next() )
    if ( !visit(edge->Node()) )
      computeVertexSolution(edge->Node(),sNode,res,visit);
}

VertexCovering::VertexCovering ( Collection * col, const deque<trail> & bag )
  :col(col),bag(bag)
{
}

deque<trail> VertexCovering::computeCovering ( bool edge, double lambda )
{
  if ( !result.empty() )
    return result;

  deque<int> suffix = getVertexSuffixTree();
  map<int,NodeInfo> sNode;
  Graph3d * grp = col->GETPOBJECT("tree",Graph3d);
  Float mWe = col->GETVALUE("sum_weight",Float);
  
  // 1 - Initialize sNode
  for ( unsigned int i = 0; i < suffix.size(); ++i )
    {
      NodeInfo infoNode;
      infoNode.dynaTab = 0;

      sNode[suffix[i]] = infoNode;
    }

  // 2 - Compute vpath sets
  computeVpath(grp,bag,suffix,sNode);

  if ( edge ) // Do edge covering using the same algorithm
    for ( unsigned int i = 0; i < suffix.size() - 1; ++i )
      sNode[suffix[i]].vpath = sNode[suffix[i]].epath;

  // 3 - Complete sNode for the next step
  for ( unsigned int i = 0; i < suffix.size(); ++i )
    {
      unsigned int size = 1 << sNode[suffix[i]].vpath.size();
      sNode[suffix[i]].dynaTab = new double [ size ];
      sNode[suffix[i]].dynaTab[0] = numeric_limits<double>::max(); // Say that is infinity
	  
      // Compute the weight of each path (can be optimized)
      for ( unsigned int x = 0; x < sNode[suffix[i]].vpath.size() ; ++x )
	{
	  double weight = 0.0;
	      
	  for ( unsigned int y = 0; y < sNode[suffix[i]].vpath[x].t.size()-1; y += 2 )
	    {
	      GEdge * edge = (*grp)[sNode[suffix[i]].vpath[x].t[y]]->Neighbours();
	      for ( ; edge != 0; edge = edge->Next() )
		if ( edge->Node() == sNode[suffix[i]].vpath[x].t[y+2] )
		  {
		    weight += edge->weight;
		    break;
		  }
	    }
	      
	  unsigned int index = 1 << x;
	  weight = -weight/mWe + lambda; // \sum_p - w(p) + \lambda = -w(P) + \lambda \card(P)
	  sNode[suffix[i]].wpath.push_back(weight);
	  sNode[suffix[i]].dynaTab[index] = weight;
	}

      // Complete the dynamic array
      for ( unsigned int x = 1; x < size; x++ )
	{
	  double weight = 0.0;
	  for ( unsigned int y = 0; y < sNode[suffix[i]].vpath.size(); ++y )
	    {
	      unsigned int index = 1 << y;
	      if (( x & index ) != 0 )
		weight += sNode[suffix[i]].wpath[y];
	    }
	  sNode[suffix[i]].dynaTab[x] = weight;
	}
    }
      
  result = pathVertexCovering(suffix,sNode);

  // Free memory
  for ( unsigned int i = 0; i < suffix.size(); ++i )
    delete[] sNode[suffix[i]].dynaTab;
 
  return result;
}

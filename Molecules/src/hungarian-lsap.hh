// =========================================================================
/* Hungarian algorithm for solving the Linear Sum Assignment
   Problem with Edition (LSAPE)
   
   authors: Sebastien Bougleux and Luc Brun
   institution: GREYC UMR 6072
                CNRS - Universite Caen Normandie - ENSICAEN
 
   reference: S. Bougleux and L. Brun
              Linear Sum Assignment Problem with Edition
              GREYC Technical Report, Oct. 2015
 
   hungarianLSAP(C,n,rho)
   Compute an assignment between two sets U and V of 
   same size n, provided a n x n non-negative cost matrix C, 
   and such that the total cost sum of the assignment is miniminal.
   The resulting assignment is provided as a bijection rho:U->V.
   Time complexity 0(n^3), Memory space 0(n^2)
   -----------------------------------------------------------
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   LSAPE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the GNU General Public License,
   and a copy of the GNU Lesser General Public License, along with 
   LSAPE. If not, see <http://www.gnu.org/licenses/>.
*/
// =========================================================================

#ifndef _HUNGARIAN_LSAP_
#define _HUNGARIAN_LSAP_

#include <limits>

// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// according to min on rows and then to min on reduced columns
// mass is the number assigned elements
// -----------------------------------------------------------
template<typename IT, class DT>
void basic_reduce_assign(const DT *C, const IT &n, IT *rho, IT *varrho, DT *u, DT *v, IT &mass)
{
  IT i = 0, j;
  DT mn, val;
  bool *U = new bool[n], *V = new bool[n];
  mass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < n; j++)
    {
      const DT &c = C[j*n+i];
      if (c < mn) mn = c;
    }
    u[i] = mn;
  }
  
  // find de min of each column
  for (j = 0; j < n; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < n; i++)
    {
      val = C[j*n+i] - u[i];
      if (val < mn) mn = val;
    }
    v[j] = mn;
    V[j] = false;
    varrho[j] = -1;
  }
  
  // assign
  for (i = 0; i < n; i++)
  {
    U[i] = false;
    rho[i] = -1;
    for (j = 0; j < n; j++)
    {
      val = C[j*n+i] - u[i] - v[j];
      if (U[i] == false && V[j] == false && val == 0)
	{ rho[i] = j; varrho[j] = i; U[i] = true; V[j] = true; mass++; break; }
    }
  }
  delete[] V; delete[] U;
}

// -----------------------------------------------------------
// Construct an alternating tree 
// and transform perform dual updates
// until an augmenting path is found
// -----------------------------------------------------------
template<typename IT, class DT>
void augment(const IT &k, const DT *C, const IT &n, const IT *rho, 
	     IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi)
{
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0;
  *SV = -1;
  zi = -1;
  
  for (i = 0; i < n; i++) { pi[i] = mx; U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      cred = C[j*n+*uluptr] - u[*uluptr] - v[j];
      if (cred < pi[*uluptr])
      {
	pred[*uluptr] = j;
	pi[*uluptr] = cred;
	if (cred == 0) { r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop; }
      }
    }
    if (lusutop == ulutop) // dual update
    {
      delta = mx;
      for (uluptr = ulutop; uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      for (uluptr = ulutop; uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == 0)  { r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop; }
      }
    } // end dual update
    i = *lusutop; ++lusutop; // i is now in SU
    if (rho[i] == -1) { zi = i; break; }
    j = rho[i];
  }
}

// -----------------------------------------------------------
// Main Hungarian algorithm
// -----------------------------------------------------------
template<typename IT, class DT>
void hungarianLSAP(const DT *C, const IT &n, IT *rho, DT *u, DT *v)
{
  IT nass = 0, i = -1, j = -1, r = -1, k = -1;
  IT *SV = NULL, *pred = NULL, *varrho = new IT[n], *U = NULL;
  DT *pi = NULL;
  
  // initialization -----------------------------------------------
  basic_reduce_assign(C,n,rho,varrho,u,v,nass);

  // augmentation of columns --------------------------------------
  if (nass < n)
  {
    U = new IT[n+1]; SV = new IT[n+1];  pi = new DT[n+1]; pred = new IT[n];
    for (k = 0; k < n; k++)
      if (varrho[k] == -1)
      {
	augment(k,C,n,rho,U,SV,pred,u,v,pi,i);  // augment always finds an augmenting path
	for (j = -1; j != k;)  // update primal solution = new partial assignment
	{
	  j = pred[i]; rho[i] = j;
	  r = varrho[j]; varrho[j] = i; i = r;
	}
      }
    delete[] U; delete[] pred; delete[] pi; delete[] SV;
  }
  delete[] varrho;
}
// -----------------------------------------------------------
#endif

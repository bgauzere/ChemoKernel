/** -----------------------------------------------------------
 * Matlab interface to the:
 * Hungarian algorithm for solving the Linear Sum Assignment Problem (LSAP)
 * authors: Sebastien Bougleux and Luc Brun
 * institution: GREYC UMR 6072
 *              CNRS - Universite Caen Normandie - ENSICAEN
 *
 * [rho,u,v] = hungarianLSAP(C)
 * Compute an assignment between two sets
 * U and V of size n, provided a nxn cost matrix C, 
 * and such that the total cost sum of the assignment is miniminal.
 *
 * The resulting assignment is provided as a mapping rho:U->V
 *
 * execute matlab file 'compile_mex' to compile this function
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
#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsap.hh>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  if (nrhs != 1 || nlhs != 3) mexErrMsgTxt("USAGE: [rho,u,v]=hungarianLSAP(C)");
  //------------------------------------------------------------------
  // 1st input : edit cost matrix
  const double *C = mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  if (nrows != ncols) mexErrMsgTxt("USAGE: cost matrix C must be square.");
  //------------------------------------------------------------------
  // outputs
  int dimr[2] = { nrows, 1 };
  plhs[0] = mxCreateNumericArray(2, dimr, mxINT32_CLASS, mxREAL);
  int *rho = (int*)mxGetData(plhs[0]);
  //------------------------------------------------------------------
  plhs[1] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
  double *u = mxGetPr(plhs[1]);
  plhs[2] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
  double *v = mxGetPr(plhs[2]);
  hungarianLSAP<int,double>(C,nrows,rho,u,v);
}


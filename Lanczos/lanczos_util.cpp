/*
 *  Copyright (C) 2010 Anthony D. Dutoi, Yasen Velkov, Michael Wormit
 *
 *  This file is part of Lanczos version 0.1.
 *
 *  Lanczos 0.1 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Author contact:  anthony.dutoi@pci.uni-heidelberg.de
 */

#include <new>
#include <complex>
#include "lanczos_util.h"

// Tony's OS X system
//#define dsbev dsbev
//#define zhbev zhbev
// Yasen's Linux system
#define dsbev dsbev_
#define zhbev zhbev_
extern "C" void dsbev(char*, char*, int*, int*,              double*,  int*, double*,              double*,  int*,              double*,           int*);
extern "C" void zhbev(char*, char*, int*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, std::complex<double>*, double*, int*);

// Declaration of Tarantelli's routines:
extern "C" void bnd2td_(int*, int*, int*, double*, double*, double*, double*, double*);
extern "C" void tddiag_(int*, int*, double*, double*, double*, int*);

namespace LanczosProjection {



float       real(const       float& x) {return x;}
double      real(const      double& x) {return x;}
long double real(const long double& x) {return x;}
float       conj(const       float& x) {return x;}
double      conj(const      double& x) {return x;}
long double conj(const long double& x) {return x;}





// This diagonalizes a real symmetric matrix which is assumed to be banded (where only the columns starting with the diagonal and 'band' 
// subdiagonals are specified in the 'dim' rows of 'mat'), and returns the eigenvalues and eigenvectors.  The eigenvectors are rows 
// (C order) of this matrix upon exit, meaning that 'mat' has dim*dim storage associated with it, but, on input, only the first
// (band+1)*dim elements are used.
//
int band_herm_diag(double *evals, double *mat, int dim, int band)
	{
	int info;

	if (dim==0) {return 0;}

	int    b1     = band + 1;
	double *bnds  = new double [b1*dim];
	double *rwork = new double [ 3*dim];

	for (int ii=0; ii<b1*dim; ii++) {bnds[ii] = mat[ii];}

	dsbev("V","L",&dim,&band,bnds,&b1,evals,mat,&dim,rwork,&info);

	delete [] rwork;  delete [] bnds;

	return info;
	}



int band_herm_diag(double *evals, std::complex<double> *mat, int dim, int band)
	{
	int info;

	if (dim==0) {return 0;}

	int                   b1    = band + 1;
	std::complex<double>* bnds  = new std::complex<double>   [b1*dim];
	double*               rwork = new              double    [ 3*dim];
	std::complex<double>* cwork = new std::complex<double>   [   dim];

	for (int ii=0; ii<b1*dim; ii++) {bnds[ii] = mat[ii];}

	zhbev("V","L",&dim,&band,bnds,&b1,evals,mat,&dim,cwork,rwork,&info);

	delete [] cwork;  delete [] rwork;  delete [] bnds;

	return info;
	}


// Diagonalizes the band-diagonal matrix using Tarantelli's routines,
// where only the the main-space components (the top part) and
// the components related to the residuals (the bottom part) of the short Lanczos vectors are obtained.
// This makes the time of diagonalizing (most probably) linearly dependent on the size of the matrix.
 int band_sym_diag_fast(double *evals, double *mat, int dim, int band)
       {
       int info;
    
       if (dim==0) {return 0;}
    
       int b1 = band + 1, nm = 2*band;
       double *e = new double[dim], *e2 = new double[dim],  *a = new double[b1*dim];
	
       for(int i = 0 ; i < b1; ++i)
	 for(int j = 0; j < dim-i; ++j) 
	   a[i+j+dim*(b1-i-1)] = mat[i+j*b1];

       // Tarantelli's routines:
       bnd2td_(&nm,&dim,&b1,a,evals,e,e2,mat);
       tddiag_(&nm,&dim,evals,e,mat, &info);

       delete [] a; delete [] e; delete [] e2;

       return info;
  }
  


}	// End namespace LanczosProjection

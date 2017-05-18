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

#ifndef _LANCZOSUTILH_
#define _LANCZOSUTILH_

#include <new>
#include <complex>

namespace LanczosProjection {



// This should make it so that real() and conj() are defined for both real and complex types inside the
// LanczosProjection namespace (or by using the function LanczosProjection::real(), etc, explicitly).
// Maybe I should be using abs in a lot of places where conj(x)*x is used ... must also be made available
// for "exotic" types
using std::real;
using std::conj;
float       real(const       float& x);
double      real(const      double& x);
long double real(const long double& x);
float       conj(const       float& x);
double      conj(const      double& x);
long double conj(const long double& x);



// The band-diagonalize wrapper necessary for Lanczos with human-readable, non-FORTRANy arguments.
int band_herm_diag(double* evals,              double*  mat, int dim, int band);
int band_herm_diag(double* evals, std::complex<double>* mat, int dim, int band);
int band_sym_diag_fast(double* evals,          double*  mat, int dim, int band);

}	// End namespace LanczosProjection

#endif	// End #ifdef _LANCZOSUTILH_

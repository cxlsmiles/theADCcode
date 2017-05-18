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

#ifndef _LANCZOSH_
#define _LANCZOSH_

#include <new>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <typeinfo>
#include "lanczos_util.h"		// Utility to make scratch memory management easier here and wrappers to do band diagonalization
#include "lanczos_engine.h"		// An essentially private member class that does the math, the stuff here does the admin
#include "lanczos_comments.h"		// Just so the compiler complains in case we ever lose track of this file, which is way too detailed for a comment

#define  RESIDUAL_DEFAULT 1e-5		// Threshold for eigenvalue-normalized residual norms for convergence
#define   OVERLAP_DEFAULT 1e-10		// Threshold for tolerated deviations in the expected (unit) overlap matrix
#define  SUBSPACE_DEFAULT 100000	// The default size of a couple of pointer arrays, indicates the maximum size of the subspace
#define            MCHUNK 10000		// The minimum "extra" amount to allocate on 'Resizeable' reallocations (affects performance/memory-use only)
#define       FREE_MEMORY true		// Just an alias for some optional arguments to make code easier to read
#define      THRESH_DUMMY 0e0		// Since this is an impossible threshold, it can be used as a signal
#define       BLOCK_DUMMY (-1)		// Since this is an impossible block size, it can be used as a signal
#define     NUMVECS_DUMMY 0             

namespace LanczosProjection {

enum {CRIT_DUMMY,CRIT_STARTPROJNORM,CRIT_EIGENVALUE};	// Different criteria by which to select those vectors that will be tested for convergence



// This Lanczos class controls the LanczosEngine class to do block Lanczos jobs ... See the stuff in lanczo_scomments.h
//
template <typename LancNum, class LancVec, class LancMat>
class Lanczos
	{

	public:

	 // Constructor/Destructor
	  Lanczos(LancMat& matrix, int criterion=CRIT_STARTPROJNORM, double threshold=RESIDUAL_DEFAULT, bool lanczosbasis=false, bool fastbanddiag = false, int subspacemax=SUBSPACE_DEFAULT);
	~Lanczos();

	 // User functions
	 void      AddStartVector(LancVec* vec);
	 void      AddStartBlock(LancVec** vecs, int nn);
	 LancNum** GetSubdiags(void);
	 LancVec** GetLanczosBasis(void);
	 int       Iterate(void);
	 int       Diagonalize(bool w_approx=false);
	 bool      DiagConverged(bool print=false, double threshold=THRESH_DUMMY, int criterion=CRIT_DUMMY);
	 double*   GetEigenvals(void);
	 LancNum*  GetEigenvecs(void);
	 double*   GetResidualNorms(void);
	 void      Lanczos2Eigen(LancNum** vecs, int nn, bool free_temp=false);
	 void      Eigen2Lanczos(LancNum** vecs, int nn, bool free_temp=false);
	 void      BuildVectors(LancVec** vecs, LancNum** coeffs, int nn, int len, int modulus=1);
	 void      BuildEigenvecs(double* eigenvals, LancVec** eigenvecs, int nn=NUMVECS_DUMMY, int criterion=CRIT_DUMMY);
	 void      Reset(LancVec **vecs, int newblock=BLOCK_DUMMY);
	 void      CleanUp(void);
	 bool      CheckOverlap(LancVec** vecs, int nn, char* title=NULL, bool printdev=true, double threshold=OVERLAP_DEFAULT);
	 bool      Orthonormalize(LancVec** vecs, int nn);

	private:

	 // These are all constant input parameters by the time 'Iterate' is called.
	 LancMat*                                mat;		// A pointer to the matrix that we are making a projection of
	 bool                                    lancbas;	// Whether we are saving the Lanczos basis
	 bool                                    fastdiag;      // Whether to speed up the diagonalization by obtaining only partial eigenvectors. Namely, just the first 'block'-sized and the last 'block'-sized part of the vectors are obtained. When this option is selected all methods requiring full eigenvectors are 'disabled'. 
	 int                                     block;		// Number of starting vectors
	 int                                     subspmax;	// Maximum number of vectors allowed in the Krylov subspace (technical point only)

	 // This is the primary information filled in by the Lanczos engine (Lanczos vectors and projection).
	 LanczosEngine<LancNum,LancVec,LancMat>* engine;	// This contains the information on the current state of the projection, and how to continue
	 LancVec**                               start;		// Pointers to starting block vectors (and generated Lanczos vectors, depending on 'lancbas')
	 LancNum**                               subdiags;	// The nonzero elements of the projected matrix are stored here
	 int                                     dim;		// Number of vectors in the projection subspace so far

	 // This keeps track of the most recent diagonalization ...
	 std::vector<double>                      evals;		// The dimd-length array of the eigenvalues of the most recently diagonalized projection
	 std::vector<LancNum>                     evecs;		// The dimd x dimd matrix of eigenvectors to the most recently diagonalized projection
	 int                                     dimd;		// Dimension of the most recently diagonalized projection, negative if approximation was used
	 int                                     diag_err;	// LAPACK-info-style code from most recent diagonalization
	 bool					 reset;		// Whether the Lanczos has been reset since the last diagonalization

	 // ... and the most recent diagonalization checked for convergence
	 double                                  thresh;	// Threshold for convergence of vector residual (argument>previous>constructor>#define)
	 int                                     crit;		// Indicates the criterion by which we will select eigenvectors to test for convergence
	 bool					 checked;	// Whether the most recent diagonalization has been checked for convergence
	 bool					 converged;	// Whether the selected eigenvectors of the most recently checked diagonalization converged
	 int*                                    list;		// List of indices of those eigenvectors selected to be tested for convergence
	 double*                                 resnorms;	// Residuals from last convergence check, for those vectors indicated in 'list'
	 int                                     oldblock;	// The block size at the last convergence check, before any potential reset

	 // Private utilities (no error checking at this level)
	 void   AddBlock(void);
	 void   MatVec(LancNum *out, LancNum *in, int nn, double eval=0);
	 double ResidualNorm(int ii, bool free_residual=false);
	 double StartProjNorm(int ii);
	 void   ListByCriterion(bool eval_sort=true, bool free_lists=false);

	  //auxiliary storage
	  std::vector<LancNum> temp;
	  std::vector<int> tempi;
	  std::vector<double> tempd;

	};





// // // // // // // // // // // // // // // // Public Functions // // // // // // // // // // // // // // // // // //



// CONSTRUCTOR
// Initialization from arguments, and zeros out other variables.  Arrays of pointers are allocated for the subdiagonals and lanczos basis.
// The only mandatory argument is the matrix object whose projection will be made.  The next two arguments specify the criterion for 
// selecting the vectors to check for convergence and the convergence threshold for the vector residual norms, for which the defaults
// are given in the above declarations.  (These can also be specified when convergence is checked.)  A boolean can be used to turn on 
// saving of the Lanczos basis, and a final integer could be given if, for some reason the default maximum subspace size is insufficient.
// (These last two arguments are considered practically debugging options).
//
template <typename LancNum, class LancVec, class LancMat>
  Lanczos<LancNum,LancVec,LancMat>::Lanczos(LancMat& matrix, int criterion, double threshold, bool lanczosbasis, bool fastbanddiag, int subspacemax) :
	mat(      &matrix      ),
	crit(     criterion    ),
	thresh(   threshold    ),
	lancbas(  lanczosbasis ), 
        fastdiag( fastbanddiag ),
	subspmax( subspacemax  )
	{
	subdiags = new LancNum* [subspmax+1];  start = new LancVec* [subspmax+1];
	for (int ii=0; ii<=subspmax; ii++) {subdiags[ii] = NULL;  start[ii] = NULL;}	// Leaves at least one NULL pointer at the end, as a signal
	evals.reserve(MCHUNK);  evecs.reserve(MCHUNK); 	// Can't do in constructor because looks like function in class def?
	temp.reserve(MCHUNK);   tempd.reserve(MCHUNK);
	block = oldblock = dim = dimd = diag_err = 0;
	converged = checked = reset = false;
	resnorms  = NULL;
	list      = NULL;
	engine    = NULL;
	if (fastdiag && (typeid(LancNum) != typeid(double))) 
	{ fastdiag = false; printf("Fast band-matrix diagonalization available for real matrices only.\n");}
	if (fastdiag && (crit == CRIT_STARTPROJNORM)) 
	{ fastdiag = false; printf("Fast band-matrix diagonalization not compatible with CRIT_STARTPROJNORM convergence criterion.\n");}
	}



// DESTRUCTOR
// This frees everything that was allocated in the course of instatiating and using this class.  It calls a function that can also be
// called externally.
//
template <typename LancNum, class LancVec, class LancMat>
Lanczos<LancNum,LancVec,LancMat>::~Lanczos() {CleanUp();}



// AddStartVector
// Takes a pointer to a vector and sets the next member of the Lanczos vector array to point to that vector.
// This function may be called many times sequentially, to add successive vectors, but it should not be called after Iterate is called,
// even if Reset has been called.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::AddStartVector(LancVec* vec) {start[block++] = vec;  return;}



// AddStartBlock
// Takes an array of pointers to vectors and sets several members of the Lancoz vector array to point to them.
// This may be called sequentially, along with AddStartVector, but it also should not be called after Iterate is called, even if Reset
// has been called.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::AddStartBlock(LancVec** vecs, int nn) {for (int ii=0; ii<nn; ii++) {AddStartVector(vecs[ii]);}  return;}



// GetSubdiags
// This returns a pointer to the raw subdiagonals array, which the user must understand the internal ordering of to use, and should
// not modify.
//
template <typename LancNum, class LancVec, class LancMat>
LancNum** Lanczos<LancNum,LancVec,LancMat>::GetSubdiags(void) {return subdiags;}



// GetLanczosBasis
// This returns a pointer to the raw start array, which will contain the Lanczos basis if these vectors are stored, otherwise it will return
// NULL.  The user should not modify the contents of this array, since some of these vectors are pointed to by other functions, including maybe
// the calling function.
//
template <typename LancNum, class LancVec, class LancMat>
LancVec** Lanczos<LancNum,LancVec,LancMat>::GetLanczosBasis(void) {return lancbas ? start : NULL;}



// Iterate
// This is the heart of this class, calling the 'LanczosEngine' class to do the actual math.
// It allocates more memory for another column of the projection (if not allocated from a previous run) and Lanczos vector (if
// saved and not already allocated) and then builds the next vector.
//
template <typename LancNum, class LancVec, class LancMat>
int Lanczos<LancNum,LancVec,LancMat>::Iterate(void)
	{
	if (block == 0) {printf("\n\nError:  Lanczos projection routine requested to iterate without a starting block.\n\n");  return -1;}
	
	if (dim+block <= subspmax)
		{if (dim == 0)		// Formalizes the start vectors as the beginning of the basis and builds what we know of the projection (ie nothing)
			{for (int ii=0; ii<block; ii++) {if (!subdiags[ii]) {subdiags[ii] = new LancNum [block+1];}}
			if (!engine) {engine = new LanczosEngine<LancNum,LancVec,LancMat> (mat, block, start, subdiags, lancbas);}}
		 else {AddBlock(); engine->AddBlock();}}
	else {printf("\n\nError:  Maximum Lanczos subspace size exceeded.\n\n");  return -1;}

	return dim+=block;
	}



// Diagonalize
// This is the most common operation one would perform on the Lanczos projection.  It diagonalizes the projection as it is up to the
// point when this is called.  One is given the option to diagonalize only the exact part or to include the approximated elements as
// well (a bit bigger matrix).  One very important consquence is the setting of the class parameter 'dimd'.  The absolute value of 'dimd'
// is the dimension of the matrix most recently diagonalized, but it is negative of the approximation was used.  This function is
// smart enough to realize if the matrix (and whether the approximation is on or off) is the same as it was the last time a diagonalization 
// was requested.  So it would do nothing the second time if called twice in a row with the same argument (but the return value and state
// of the object will be as though it had).  It also sets checked to false, since this new diagonalzation has not been checked for convergence
// and it set reset to false because the eigenvectors/values do belong to some submatrix of the projection again.
//
template <typename LancNum, class LancVec, class LancMat>
int Lanczos<LancNum,LancVec,LancMat>::Diagonalize(bool w_approx)
	{
	// First figure out what the dimension of the diagonalization is, based on dim and whether we us an approximation (denoted by negative dimd)
	int dimd_copy = dimd;
	dimd = w_approx ? -dim : std::max(dim-block,0);
	oldblock = block;
	if (!fastdiag) {
	// If the exact same dimension (and yes or no on approximation), do nothing.  If a reset was done, we should diagonalize, because dimd==dimd_copy
	// is clearly a coincidence from having called Iterate just the right number of times.
	if (dimd!=dimd_copy || reset)
		{evals.resize(abs(dimd));  evecs.resize((unsigned)std::max(abs(dimd),block+1)*abs(dimd));
		 for (unsigned int ii=0; ii<abs(dimd); ii++) {for (unsigned int jj=0; jj<=block; jj++) {evecs[ii*(block+1)+jj] = subdiags[ii][jj];}}
		 diag_err = band_herm_diag((double *)&evals.front(), (LancNum *)&evecs.front(), abs(dimd), block);
		 checked = reset = false;}
	} else {//fastdiag:                the size of the storage of the partial vectors is dimd*2block, big enough to fit subdiags
	        {evals.resize(abs(dimd));  evecs.resize((unsigned)std::max(abs(dimd),block+1)*2*block);
		 for (unsigned int ii=0; ii<abs(dimd); ii++) {for (unsigned int jj=0; jj<=block; jj++) {evecs[ii*(block+1)+jj] = subdiags[ii][jj];}}
		 diag_err = band_sym_diag_fast((double *)&evals.front(), (double *)&evecs.front(), abs(dimd), block);
		 checked = reset = false;}
	  
	}
	return diag_err;
	}



// DiagConverged
// This tests the eigenvectors of the most recent *diagonalization* for convergence (whether or not anything has been called since then).
// Convergence is tested by comparing the residuals of selected vectors to a threshold.  The routine is smart enough to know if the most 
// recent diagonalization has been already been checked for convergence for the same selection of vectors for the same threshold, 
// and simply returns that result (and the default for both the threshold and selection criterion is set to the same as the previous
// call.)  The vectors are selected according to a user-chosen criterion, 'crit', the most common being to select 'block' vectors 
// comprising the lowest eigenvectors or those most similar to the starting block.  An important caveat is that the current projection 
// must be of a dimension at least 'block' larger than the vectors being tested to compute residuals.  Therefore, if the diagonalization 
// used an approximate matrix, 'Iterate' must be called again before the convergence may be checked.  One may access information from the
// previous diagonalization even if Reset has been called, so long as different vectors or threshold are not requested, and that 
// diagonalization had been checked before the call to Reset.  This routine always refers to the block size as 'oldblock', because, if
// Reset has been called and no new diagonalization has been done, all operations should refer to the block size at the time of the
// diagonalization that preceded the reset.  Calling the diagonalizer sets the variable oldblock to the value of block.
//
template <typename LancNum, class LancVec, class LancMat>
bool Lanczos<LancNum,LancVec,LancMat>::DiagConverged(bool print, double threshold, int criterion)
	{

	  if (threshold == THRESH_DUMMY) {threshold = thresh;}		// Wanted to have default arguments like criterion=crit
	if (criterion ==   CRIT_DUMMY) {criterion = crit;}		// above, but apparently this is illegal!

	if (print) {printf("diag_dim: %4d  residual_thresh: %.1le     ",abs(dimd),threshold);}

	// In this case there a not even enough eigenvectors to check for convergence ... most likely situation is that this was called
	// after the first round of Iterate and Diagonalize, where there was nothing to do, and nothing to check.
	if (abs(dimd) < oldblock)
		{if (print) {printf("Diagonalization too small for convergence check.\n");  fflush(stdout);}	
		 checked = true;  thresh = threshold;  crit = criterion;
		 return converged=false;}

	// Do the convergence checking.  Three things could cause us to actually check, rather than just returning the previous answer.  
	// Either we have not checked yet since the last diagonalization, or we have checked but we are crurious about different vectors 
	// or a tighter/looser threshold.  If only the threshold is new, we can just use saved information, but if either of the other 
	// conditions are true we need new residuals, so we check that there is enough matrix projection to compute them and that the 
	// eigenvectors are from that projection.
	double max=-1, avg_eval=0;
	if (!checked || thresh!=threshold || crit!=criterion)
		{if ((!checked || crit!=criterion) && (reset || dim-abs(dimd)<oldblock))
			{printf("\n\nError:  Cannot compute residuals.  Lanczos was reset or projection is too small.\n\n");  fflush(stdout);
			 return converged=false;}
		 else if (!checked || crit!=criterion)
			{checked = true;  crit = criterion;
			 ListByCriterion();  if (!resnorms) {resnorms = new double [oldblock];}
			 for (int ii=0; ii<oldblock; ii++) {resnorms[ii] = ResidualNorm(list[ii]);}}
		 for (int ii=0; ii<oldblock; ii++) {if (resnorms[ii]>max) {max = resnorms[ii];}}
		 for (int ii=0; ii<oldblock; ii++) {avg_eval += fabs(evals[list[ii]]);}  avg_eval /= oldblock;
		 thresh = threshold;
		 converged = max/avg_eval < thresh;}

	// Print stuff (optional)
	if (print)
		{if (max<0)
			{for (int ii=0; ii<oldblock; ii++) {if (resnorms[ii]>max) {max = resnorms[ii];}}
			 for (int ii=0; ii<oldblock; ii++) {avg_eval += fabs(evals[list[ii]]);}  avg_eval /= oldblock;}
		
		 printf("max_residual/avg_eigenvalue: %.1le  ", max/avg_eval);
		 if (!converged)
			{printf("unconverged\n");  
//for (int ii=0; ii<oldblock; ii++) {printf("%.2le  ",evals[list[ii]]);}  printf("\n");
//for (int ii=0; ii<oldblock; ii++) {printf("%.2le  ",resnorms[ii]);}     printf("\n\n");
			 fflush(stdout);}
		 else
			{printf("  converged!\n");
//for (int ii=0; ii<oldblock; ii++) {printf("%.2le  ",evals[list[ii]]);}  printf("\n");
			 printf("Raw residual norms (average eigenvalue = %.2le):\n",avg_eval);
			 for (int ii=0; ii<oldblock; ii++) {printf("%.2le  ",resnorms[ii]);}  printf("\n\n");  fflush(stdout);}}

	return converged;
	}



// GetEigenvals
// Returns a pointer to the array of eigenvalues from the most recent diagonalization, which the user should not modify the contents of.
// The user must know what dimd is to use this, but that isn't hard, since Iterate returns dim.
//
template <typename LancNum, class LancVec, class LancMat>
double* Lanczos<LancNum,LancVec,LancMat>::GetEigenvals(void) {return (double *)&evals.front();}



// TODO: Call this method in DiagConverged
// Revision 2.0 Mon Dec  7 17:38:37 CET 2009  yasen
// GetResidualNorms()
// Returns a pointer to the array of residual norms from the most recent diagonalization.
//
template <typename LancNum, class LancVec, class LancMat>
double*  Lanczos<LancNum,LancVec,LancMat>::GetResidualNorms(void)
        {
	  if (resnorms) delete [] resnorms;
	  
	  ListByCriterion();
	  resnorms = new double [abs(dimd)];

	  for (int ii=0; ii<abs(dimd); ii++) {resnorms[ii] = ResidualNorm(list[ii]);}
	  
	  return resnorms;
	}



// Lanczos2Eigen
// The user must know what dimd is to use this, but that isn't hard, since Iterate returns dim.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::Lanczos2Eigen(LancNum** vecs, int nn, bool free_temp)
	{
	  //static std::vector<LancNum> temp; temp.reserve(MCHUNK);

	if (free_temp) {temp.clear();  return;}
	if (fastdiag)  {printf("\n Can't run Lanczos2Eigen with partial eigenvectors\n"); return;}

	temp.resize(abs(dimd));
	for (int ii=0; ii<nn; ii++)
		{for (int jj=0; jj<abs(dimd); jj++) {temp[jj] = 0;  for (int kk=0; kk<abs(dimd); kk++) {temp[jj] += vecs[ii][kk]*evecs[jj*abs(dimd)+kk];}}
		 for (int jj=0; jj<abs(dimd); jj++) {vecs[ii][jj] = temp[jj];}}

	return;
	}



// Eigen2Lanczos
// The user must know what dimd is to use this, but that isn't hard, since Iterate returns dim.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::Eigen2Lanczos(LancNum** vecs, int nn, bool free_temp)
	{
	  //static std::vector<LancNum> temp;temp.reserve(MCHUNK);

	if (free_temp) {temp.clear();  return;}
	if (fastdiag)  {printf("\n Can't run Eigen2Lanczos with partial eigenvectors\n"); return;}

	temp.resize(abs(dimd));
	for (int ii=0; ii<nn; ii++)
		{for (int jj=0; jj<abs(dimd); jj++) {temp[jj] = 0;}
		 for (int kk=0; kk<abs(dimd); kk++) {for (int jj=0; jj<abs(dimd); jj++) {temp[jj] += vecs[ii][kk]*evecs[kk*abs(dimd)+jj];}}
		 for (int jj=0; jj<abs(dimd); jj++) {vecs[ii][jj] = temp[jj];}}

	return;
	}



// New func!
// The user must know what dimd is to use this, but that isn't hard, since Iterate returns dim.
//
template <typename LancNum, class LancVec, class LancMat>
LancNum* Lanczos<LancNum,LancVec,LancMat>::GetEigenvecs(void)
	{
	  return (LancNum *) &evecs.front();
	}



// BuildVectors
// This builds nn vectors, whose len Lanczos-basis coefficients are given in coeffs, in the representation 
// of the original matrix for which the projection is being computed.  This leaves the object in such a state 
// as it was before, if the Lanczos basis was saved, or as if it had been Reset (see below) and then had 
// Iterate called on it len/block times, without ever touching the Diagonalizer.
// The final optional argument is passed all the way through to the matrix-vector multiplication, allowing for
// the more efficient generation of the lanczos basis by acting on more than one vector at a time.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::BuildVectors(LancVec** vecs, LancNum** coeffs, int nn, int len, int modulus)
	{
	// Check that we have either called reset or stored the Lanczos basis.
	if (dim!=0 && !lancbas) {for (int ii=0; ii<nn; ii++) {vecs[ii] = NULL;}  return;}

	// Pass it off to the engine (which assumes that we have stored the basis or called a reset).
	
	engine->BuildVecs(vecs,coeffs,nn,len,modulus);
		
	// Leave the object in a self-consistent state, and let it know that the projection and eigenvectors may not match.
	if (!lancbas) {dim = len;  reset = true;}

	return;
	}


// Revision 1.0 Wed Nov 18 14:09:05 CET 2009  yasen
// BuildEigenvecs generates a number of vectors selected by the user: block replaced with nn
//
// BuildEigenvecs
// This builds the eigenvectors of the most recent diagonalization, which are selected by the criterion mentioned in DiagConverged.  It 
// builds these in the representation of the original matrix for which the projection is being computed.  This leaves the object in such
// a state as it was before, if the Lanczos basis was saved, or as if it had been Reset (see below) and then had Iterate called on it
// abs(dimd)/block times, without ever touching the Diagonalizer.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::BuildEigenvecs(double* eigenvals, LancVec** eigenvecs, int nn, int criterion)
	{
	if (fastdiag)  {printf("\n Can't run BuildEigenvecs with partial eigenvectors\n"); return;}

	if (criterion == CRIT_DUMMY) {criterion = crit;}		// See comment near similar lines in DiagConverged

	if (nn == NUMVECS_DUMMY) {nn = block;}		// If no size is specified use the block size as a default

	// Check that we have at least block eigenvectors to build
	if (abs(dimd) < nn) {for (int ii=0; ii<nn; ii++) {eigenvals[ii] = NAN;  eigenvecs[ii] = NULL;}  return;}
	
	// Run this if diagonalized since last convergence check, or want different vectors
	if (!checked || crit!=criterion) {crit = criterion;  ListByCriterion();}

	// Collect the coefficients and pass it off to BuildVectors (which assumes that we have stored the basis or called a reset).
	LancNum** Eigenvecs = new LancNum* [nn];
	for (int ii=0; ii<nn; ii++) {eigenvals[ii] = evals[list[ii]];  Eigenvecs[ii] = evecs + list[ii]*abs(dimd);}
	
	BuildVectors(eigenvecs,Eigenvecs,nn,abs(dimd));
	
	delete [] Eigenvecs;

	return;
	}



// Reset
// This does a limited resetting of the object, but does not leave it in the "fresh from the constructor" state, because no allocations 
// are freed, since many allocations depend on the block size.  The object can only be used for Lanczos iterations on the same type of 
// vector class and for the same block size, but after that, the reset is rather general.  The argument specfies a new set of 'block' 
// starting vectors and Iterate may be called.  Living up to the letter of the descriptions of the functions in this class, the 
// diagonalizer information is not cleared out, and even after Reset is called the previous eigenvectors and values are available and all 
// information that is supposed to depend on the last diagonalization does exactly that, regardless of the fact that we now have a 
// smaller matrix, or even a different one.  In spite of this, one may request a new diagonalization of the new matrix, and this will 
// behave as expected.  The most important use of this function will be to reset the starting vectors before and/or after calling 
// BuildEigenVectors when the Lanczos basis is not stored, which is the primary reason for not wiping out the diagonalizer information.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::Reset(LancVec **vecs, int newblock)
	{
	if (newblock == BLOCK_DUMMY) {newblock = block;}		// See comment near similar lines in DiagConverged

	if (newblock != block)
		{if (newblock > block)
			{if (subdiags)
				{int ii = 0;  while (subdiags[ii]) {delete [] subdiags[ii];  subdiags[ii++] = NULL;}
				 for (ii=0; ii<newblock; ii++) {subdiags[ii] = new LancNum [newblock+1];}}
			 if (start)
				{int ii = block;  while (start[ii]) {ii++;}
				 for (int jj=block; jj<newblock; jj++, ii++)
					{if (ii < subspmax) {start[ii] = start[jj];}
					 else               {delete      start[jj];}
					 start[jj] = NULL;}}
			 if (resnorms)
				{double* temp = resnorms;  resnorms = new double [newblock];
				 for (int ii=0; ii<block; ii++) {resnorms[ii] = temp[ii];}
				 delete [] temp;}
			 if (list)
				{int*    temp = list;      list     = new int    [newblock];
				 for (int ii=0; ii<block; ii++) {list[ii]     = temp[ii];}
				 delete [] temp;}}
		 else
			{if (start)
				{for (int ii=newblock; ii<block; ii++) {start[ii] = NULL;}
				 int ii = block;  int jj = newblock;  while (start[ii]) {start[jj++] = start[ii]; start[ii++] = NULL;}}}
		 oldblock = block;  block = newblock;}

	for (int ii=0; ii<block; ii++) {start[ii] = vecs[ii];}		// Must do before restarting engine, because engine remembers start
	engine->Reset(block);						// No harm done if block size not changed
	reset = true;			// Signal that the previous diagonalization is no longer from a submatrix of the current projection
	dim   = 0;			// Start over!

	return;
	}



// CleanUp
// This can be called externally, and does everything the destructor should do, with the purpose of freeing up memory, the object is 
// rendered useless afterwards.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::CleanUp(void)
	{
	int ii;
	if (start   ) {ii = block;  while (   start[ii]) {delete    start[   ii++];}  delete [] start;     start    = NULL;}
	if (subdiags) {ii = 0;      while (subdiags[ii]) {delete [] subdiags[ii++];}  delete [] subdiags;  subdiags = NULL;}
	if (engine)   {delete    engine;    engine   = NULL;}
	if (resnorms) {delete [] resnorms;  resnorms = NULL;}
	if (list)     {delete [] list;      list     = NULL;}
	ResidualNorm(       0,FREE_MEMORY);
	ListByCriterion(false,FREE_MEMORY);
	Eigen2Lanczos(NULL, 0,FREE_MEMORY);
	Lanczos2Eigen(NULL, 0,FREE_MEMORY);
	evals.clear();  evecs.clear();
	block = oldblock = dim = dimd = diag_err = subspmax = 0;
	converged = checked = reset = lancbas = false;
	thresh = THRESH_DUMMY;
	crit   =   CRIT_DUMMY;
	mat = NULL;
	return;
	}



// CheckOverlap
// This checks that the set of vectors given as an argument are orthonormal, given the threshold in the final argument.
// If 'title' is not NULL, then the overlap matrix is printed with the given title (deviation from the unit matrix).
// This function depends on no class parameters, except the template types.  Therefore it need not be a member function but it is 
// convenient utility, so we keep it here.
//
template <typename LancNum, class LancVec, class LancMat>
bool Lanczos<LancNum,LancVec,LancMat>::CheckOverlap(LancVec** vecs, int nn, char* title, bool printdev, double threshold)
	{
	LancNum one;
	bool    ok;

	if (title)
		{printf("\nOverlap matrix for %svectors",title);
		 printf("%s:\n", printdev ? " (deviation from unit matrix, per-element norms of differences)" : "");}

	one = 1;
	ok  = true;
	for (int ii=0; ii<nn; ii++)
		{if (title) {printf("%*s",ii*(printdev?10:24),"");}
		 for (int jj=ii; jj<nn && (ok||title); jj++)
			{LancNum olap = (*vecs[ii])*(*vecs[jj]);
			 if (title && !printdev) {printf(" % .20lf",olap); fflush(stdout);}	// This line breaks if LancNum!=double (just comment it out)
			 if (ii==jj) {olap -= one;}
			 double norm = sqrt(real(conj(olap)*olap));
			 if (title &&  printdev) {printf("  %.2le", norm); fflush(stdout);}
			 if (norm > threshold) {ok = false;} }
		 if (title) {printf("\n");} }
	if (title) {printf("\n");}  fflush(stdout);

	return ok;
	}



// Orthonormalize
// This orthonormalizes the set of vectors given as an argument using Graham-Schmidt (later, several options will be implemented).
// It returns a quality control variable, based on how nonorthogonal the vectors were, if it is false, one might want to repeat the procedure.
// This function depends on no class parameters, except the template types.  Therefore it need not be a member function but it is 
// convenient utility, so we keep it here.
//
template <typename LancNum, class LancVec, class LancMat>
bool Lanczos<LancNum,LancVec,LancMat>::Orthonormalize(LancVec** vecs, int nn)
	{
	bool ok = true;
//printf("orthonormalize called\n");  fflush(stdout);
	for (int ii=0; ii<nn; ii++)
		{for (int jj=0; jj<=ii; jj++)
			{LancNum olap = (*vecs[jj])*(*vecs[ii]);
			 if (jj < ii) { vecs[ii]->LinCombo(1,*vecs[ii],-olap,*vecs[jj]);}
			 else         {*vecs[ii] /= sqrt(real(olap));  if (sqrt(real(olap)) < (double)1/2) {ok = false;}}}}
	return ok;
	}





// // // // // // // // // // // // // // // // Private Functions // // // // // // // // // // // // // // // // // //



// AddBlock
// If another block of vectors are added, this allocates memory for them, and the relevant subdiagonals.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::AddBlock(void)
	{
	for (int ii=0; ii<block; ii++) {if (!subdiags[dim+ii]) {subdiags[dim+ii] = new LancNum [block+1];}}
	if (lancbas) {for (int ii=0; ii<block; ii++) {if (!start[dim+ii]) {start[dim+ii] = new LancVec (**start);}}}
	return;
	}



// MatVec
// This acts the matrix projection, as stored in subdiags, on the vector 'in' of length 'nn', which must be shorter than dim by at
// least block, and stores the result to the vector 'out' of length 'nn'+block.  The entries of 'in' beyond 'nn' are simply assumed to be
// zero.  We note that, under this collection of conditions, the matrix*vector product is exact.  The final optional argument (default is 
// zero) is the value of a constant (diagonal matrix) to be subtracted from the the projection; this is useful for computing residuals of 
// approximate eigenvectors in one step, which is why the parameter bears the name 'eval'.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::MatVec(LancNum *out, LancNum *in, int nn, double eval)
	{
	LancNum Eval;  Eval = eval;

	for (int ii=0; ii<nn+block; ii++) {out[ii] = 0;}
	for (int ii=0; ii<nn;       ii++)
	  {out[ii] += (subdiags[ii][0]-Eval)*in[ii];
		 for (int jj=1; jj<=block; jj++)
		   {                    out[ii+jj] +=      subdiags[ii   ][jj] *in[ii];
			 if (ii-jj>=0) {out[ii-jj] += conj(subdiags[ii-jj][jj])*in[ii];}}}

	return;
	}



// ResidualNorm 
// This computes the norm of the residual of the 'ii'-th eigenvector of the most recent diagonalization, assuming that the dimension of
// the projection is at least block larger than the dimension of the most recently diagonalized projection.  It allocates some scratch 
// memory which can be freed by calling this with the second argument set to true.
//
template <typename LancNum, class LancVec, class LancMat>
double Lanczos<LancNum,LancVec,LancMat>::ResidualNorm(int ii, bool free_residual)
	{
	if (free_residual) {temp.clear();  return 0;}

	if (!fastdiag) {

	  temp.resize(dim);
	  MatVec((LancNum *)&temp.front(),&evecs.front()+ii*abs(dimd),abs(dimd),evals[ii]);
	
	} else {
	  // Performs just the multiplication between the band-diagonal matrix and the lowest part of the partial eigenvector
	  temp.resize(block);
	  int nn = abs(dimd);  
	  LancNum* out = &temp.front();  
	  LancNum* in = &evecs.front()+(2*ii+1)*block; //This line assumes the partial vectors are being used here.
	  
	  for (int ii=0; ii<block; ii++) {out[ii] = 0;}
	  for (int ii=nn-block; ii<nn;       ii++) {
	    int ind = ii - nn + block;
	    
	    for (int jj=0; jj<block; jj++) {
	      if (ind+jj>=block) break;
	      out[ind] +=  subdiags[ii+jj][block-jj]*in[ind+jj];
	    }
	    
	  }
	}
	
	double norm2 = 0;  for (int jj=0; jj<dim; jj++) {norm2 += real(conj(temp[jj])*temp[jj]);}

	return sqrt(norm2);
	}



// StartProjNorm
// This computes the projection of the 'ii'-th eigenvector of the most recent diagonalization into the subspace spanned by the starting block.
//
template <typename LancNum, class LancVec, class LancMat>
double Lanczos<LancNum,LancVec,LancMat>::StartProjNorm(int ii)
	{
	if (fastdiag)  {printf("\n Can't run StartProjNorm with partial eigenvectors\n"); return 0.;}

	LancNum *vec = &evecs.front()+ii*abs(dimd);
	double norm2 = 0;  for (int jj=0; jj<block && jj<abs(dimd); jj++) {norm2 += real(conj(vec[jj])*vec[jj]);}
	return sqrt(norm2);
	}



// Revision 1.0 Thu Nov 19 11:35:02 CET 2009  yasen
// ListByCriterion creates the full list of ordered (by a criterion) vectors: block replaced with abs(dimd)
//
// ListByCriterion
// This allocates (if necessary) and fills the values of the array 'list'.  It reads the parameter 'crit' which indicates a user chosen
// criterion for selecting eigenvectors of the most recent diagonalization.  If the value of 'crit' is positive, then it selects 'block' vectors 
// which have the largest value of that criterion measure (assumed to take a real value), if it is negative, it chooses the lowest.  The 
// indices of the selected vectors are the integer values that fill 'list'.  If 'eval_sort' is true, it then arranges these indices in ascending 
// order of the eigenvalues associated with the selected eigenvectors.  It should be clear that the dimension of the most recent diagonalization
// needs to be at least 'block'.
//
template <typename LancNum, class LancVec, class LancMat>
void Lanczos<LancNum,LancVec,LancMat>::ListByCriterion(bool eval_sort, bool free_lists)
	{
	  //static std::vector<int>    tempi;		// Statics retain their values between objects of the same class, so resizeable is safest.
	  //static std::vector<double> tempd;tempd.reserve(MCHUNK);

	if (free_lists) {tempi.clear();  tempd.clear();  return;}

	if (!list) {list = new int [abs(dimd)];}  tempi.resize(abs(dimd));

	if      (crit == -CRIT_EIGENVALUE) {for (int ii=0; ii<abs(dimd); ii++) {tempi[ii] =             ii;}}	// Treat this special case most efficiently
	else if (crit ==  CRIT_EIGENVALUE) {for (int ii=0; ii<abs(dimd); ii++) {tempi[ii] = abs(dimd)-1-ii;}}	// Treat this special case most efficiently
	else
		{tempd.resize(abs(dimd));
		 if (abs(crit) == CRIT_STARTPROJNORM) {for (int jj=0; jj<abs(dimd); jj++) {tempd[jj] = (crit<0?-1:1)*StartProjNorm(jj);}}
		 // other else-if statements can be added here, or maybe switch-break is better?
		 for (int ii=0; ii<abs(dimd); ii++)
			{double max       = std::numeric_limits<double>::min()/2.; 	// = -HUGE_VAL/2;
			 for (int jj=0; jj<abs(dimd); jj++) {if (tempd[jj] > max) {max = tempd[jj];  tempi[ii] = jj;}}
			 tempd[tempi[ii]] = std::numeric_limits<double>::min();}}	// = -HUGE_VAL;

	if (eval_sort && crit!=-CRIT_EIGENVALUE)	// Handles the most common special case most efficiently
		{tempd.resize(abs(dimd));
		 for (int ii=0; ii<abs(dimd); ii++) {tempd[ii] = evals[tempi[ii]];}
		 for (int ii=0; ii<abs(dimd); ii++)
			{int kk;
			 double min = std::numeric_limits<double>::max()/2.;	// = HUGE_VAL/2;
			 for (int jj=0; jj<abs(dimd); jj++) {if (tempd[jj]<min) {min = tempd[jj];  list[ii] = tempi[jj];  kk = jj;}}
			 tempd[kk]  = std::numeric_limits<double>::max();}	// = HUGE_VAL;
		 bool ok = true;  for (int ii=0; ii<abs(dimd) && ok; ii++) {if (list[ii] >= abs(dimd)) {ok = false;}}
		 if (!ok) {printf("\n\nwarning:  There are lower eigenvectors outside of the set tested for convergence.\n\n");  fflush(stdout);}}
	else {for (int ii=0; ii<abs(dimd); ii++) {list[ii] = tempi[ii];}}

	return;
	}


}	// End namespace LanczosProjection




#endif	// End #ifdef _LANCZOSH_

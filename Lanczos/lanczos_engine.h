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

#ifndef _LANCZOSENGINEH_
#define _LANCZOSENGINEH_

#include <new>
#include <algorithm>
#include <cmath>
#include "lanczos_util.h"	// Defines conj() and real() for both real and complex types

namespace LanczosProjection {



// This class is the private workhorse on which the Lanczos class depends.
//
template <typename LancNum, class LancVec, class LancMat>
class LanczosEngine
	{

	public:

	 LanczosEngine(LancMat *matrix, int blocksize, LancVec **startvecs, LancNum **subdiagcols, bool lanczosbasis);
	~LanczosEngine();

	 void Reset(int newblock=0);
	 void AddVectors(int nn);
	 void AddVector();
	 void AddBlock();
	 void BuildVecs(LancVec **toBuild, LancNum **coefficients, int numvecs, int lenvecs, int modulus=1);

	private:

	 LancMat*        mat;			// The almighty matrix, acts on the vectors to generate the Lanczos projection of itself
	 int             block;			// The number of starting vectors
	 LancVec**       start;			// Pointers to the starting block vectors (and Lanczos vectors, if stored)
	 LancNum**       subdiags;		// The elements of the projected matrix will be stored here
	 int             dim;			// The dimension of exact part of the projection so far
	 LancVec**       vectors;		// Manage pointers to vectors currently being worked on
	 LancVec**       scratch;		// Pointers to scratch space (if needed, ie, not storing Lanczos vectors)

	 void Initialize(bool lanczosbasis);
	 void SubdiagClear(int ii);
	 void ShiftPointers(int nn);

	};





// // // // // // // // // // // // // // // // Public Funtions // // // // // // // // // // // // // // // // // //



// CONSTRUCTOR
// The data section sets those things that are constants during the iterations that build the projected matrix.
// The first argument 'matrix' is an object that somehow represents the operator we are making a projection of
// via its action on vectors.  An array of pointers to starting vectors is given in 'startvecs', which should
// point to at least 'blocksize' vectors, and exactly 'blocksize' number of vectors should contain meaningful
// information, being a set of orthonormal vectors (the remaining vectors may recieve output).  If 'lanczosbasis'
// is true, the Lanczos basis is to be stored, and 'startvecs' should contain additional pointers to the locations
// where the Lanczos vectors being generated are to be stored (or these pointers will be assigned by the calling
// function by the time that Lanczos vector is to be stored).  The starting vectors are the first 'blocksize' Lanczos 
// vectors by definition.  An array of pointers to output the columns of the projected matrix is given in 'subdiagcols'.  
//
// After construction, the public function 'AddVectors' is called to generate the desired number of successive columns 
// of the projection and Lanczos vectors.  The user may concurrently check the projection according to their convergence 
// criterion and add vectors to the space as desired, allocating more projection columns and space for Lanczos vectors 
// (if saved), and assigning the pointers to these allocations to the successive elements of the 'subdiagcols' and 
// 'startvecs' arrays, respectively.
//
template <typename LancNum, class LancVec, class LancMat>
LanczosEngine<LancNum,LancVec,LancMat>::LanczosEngine(LancMat *matrix, int blocksize, LancVec **startvecs, LancNum **subdiagcols, bool lanczosbasis) :
 mat(     matrix     ),
 block(   blocksize  ),
 start(   startvecs  ),
 subdiags(subdiagcols)
	{
	if (lanczosbasis)
		{scratch = start + block;}		// The start array has pointers to the output location for Lanczos vectors
	else
		{scratch = new LancVec* [block];	// If we are not saving Lanczos vectors we need scratch vectors 
		 for (int ii=0; ii<block; ii++) {scratch[ii] = new LancVec(**start);}}
	vectors = NULL;
	Initialize(lanczosbasis);
	}



// DESTRUCTOR
// This simply frees any memory that was allocated for either scratch vector space or to keep track of pointers.
//
template <typename LancNum, class LancVec, class LancMat>
LanczosEngine<LancNum,LancVec,LancMat>::~LanczosEngine()
	{
	if (scratch) {for (int ii=block-1; ii>=0; ii--) {delete scratch[ii];}  delete [] scratch;}
	delete [] vectors;
	}



// Reset
// This sets the internal state of the engine as if it had just been constructed.  The main advantage of this over
// destroying and creating a new object is that some of the setup from the previous construction is left in tact.
// For example, the pointers to 'start', 'subdiags' and 'mat' are left as is.  Also, for problems of the same block
// size, the scratch vectors are simply reused, rather than deallocated and reallocated, and, even when the block
// size is changed, some or all of the scratch vectors are reused (the restart can only handle vectors of the same
// type/class, since that is a template parameter).  While the option whether or not to store the Lanczos basis may 
// not currently be changed upon a restart (fix this, should be easy), the block size can be.  An important caveat
// is that this changes the expected lengths of the arrays pointed to by the elements of 'subdiags'.  Since 'subdiags'
// was an input parameter, the user presumably had a memory of this pointer and can make the necessary alterations
// to the content it points to.  Also, the philosophy of this engine is that *in principle* all of the memory pointed
// to by the elements of 'start' or 'subdiags' was allocated before the constructor was called, and this remains
// true of a call to this 'Reset'.  This is most important as concerns the first 'block' members of the 'start'
// and 'subdiag' arrays, as these pointers are copied into a private bookkeeping array or have their contents set
// to zero here, respectively.  While it is not exactly necessary that the information in the vectors represents 
// an orthonormal set at the time of calling this, these vectors had better be the vectors that will (or already) 
// hold the starting vectors.  As has been said before *in practice* the pointers to the later Lanczos vectors 
// (if saved) and subdiagonals may be allocated just before 'AddVectors' is expected to fill them, since they will 
// not have been referenced before this.  The user further needs to take care that the vectors in 'start' were 
// overwitten since the time of beginning the first Lanczos proecedure, unless the Lanczos basis was saved.  One 
// very efficient use of this reset is before the calling BuildVecs(), when the Lanczos basis has not been saved, 
// because then one knows that all the space for 'subdiags' that is needed has already been allocated, and one 
// just needs to remake the starting vectors in their original location before calling BuildVecs.  If the only
// argument to this function is nonzero (zero is the default), then the block size is changed to that size.
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::Reset(int newblock)
	{
	if (newblock!=0 && newblock!=block)
		{if (scratch)
			{LancVec** temp = scratch;  int ii;
			 scratch = new LancVec* [newblock];
			 for (ii=0;  ii<std::min(block,newblock); ii++) {scratch[ii] = temp[ii];}
			 for (ii=ii; ii<std::max(block,newblock); ii++)
				{if (block>newblock) {delete temp[ii];}
				 else                {scratch[ii] = new LancVec(**start);}}
			 delete [] temp;}
		 block = newblock;}

	bool lanczosbasis = !scratch;
	if (lanczosbasis) {scratch = start + block;}
	Initialize(lanczosbasis);
	return;
	}



// AddVector
// This is simply a wrapper for AddVectors
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::AddVector() {AddVectors(1);      return;}



// AddBlock
// This is simply a wrapper for AddVectors
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::AddBlock()  {AddVectors(block);  return;}



// AddVectors
// This is main event which is iterated to build the next Lanczos vectors and perform the projection of the matrix 
// into the space spanned by the vectors thus far.  It takes as an input parameter (nn<=block, not checked!!) the 
// number of new Lanczos vectors that are be added to the Lanczos set.
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::AddVectors(int nn) 
{

  LancVec**     old_set = vectors;
  LancVec** current_set = vectors +   block;
  LancVec**     new_set = vectors + 2*block;		// Logically distinct, but may be physically the same as the old vectors

  // Step 1:  Orthogonalize against the vectors in the old block (addresses in old and new blocks might be the same)
  for (int new_vec = 0;  new_vec < nn;  ++new_vec)
    
    for (int old_vec = new_vec;  old_vec < block;  ++old_vec) {

      int col = dim - block + old_vec;
      int row = dim + new_vec - col;

      if (old_set[old_vec])
         new_set[new_vec]->LinCombo(old_vec!=new_vec,*new_set[new_vec], -conj(subdiags[col][row]),*old_set[old_vec]);
      else
        *new_set[new_vec] = 0;

    }
 

  // Step 2:  The soul.  Act the matrix on the next group of vectors.
  int ii = 0;
  while (ii < nn)  ii += (*mat)(new_set+ii,current_set+ii,nn-ii);


  // Step 3:  Orthogonalize the new vectors against the current vectors and the preceding new vectors
  for (int new_vec = 0;  new_vec < nn;  ++new_vec) {

    // Step 3.1:  Orthogonalize against the "ealier" current vectors ... the subdiags elements will exist by the time they are used
    for (int current_vec = 0;  current_vec < new_vec;  ++current_vec) {

      int col = dim + current_vec;
      int row = dim + new_vec - col;

      new_set[new_vec]->LinCombo(1,*new_set[new_vec],-conj(subdiags[col][row]),*current_set[current_vec]);
      
    }

    // Step 3.2: Othogonalize against the "later" current vectors and preceding new vectors in one loop (iterates across current/new boundary)
    for (int current_vec = new_vec;  current_vec < new_vec + block;  ++current_vec) {

      int col = dim + new_vec;
      int row = dim + current_vec - col;

      subdiags[col][row] = (*current_set[current_vec])*(*new_set[new_vec]);		// These overlaps are also matrix elements (Krylov def)

      new_set[new_vec]->LinCombo(1,*new_set[new_vec], -subdiags[col][row],*current_set[current_vec]);

    }

    // Step 3.3: Generate the final matrix element and normalize the vector
    subdiags[dim+new_vec][block] = sqrt((*new_set[new_vec])||2);			// This element is the norm of what is left
    *new_set[new_vec] /= subdiags[dim+new_vec][block];

  }


  // Step 4:  Increment the dimension, do pointer admin, and fill in approximation for unknown diagonals (etc)
  ShiftPointers(nn);
  for (ii = dim;  ii < dim+nn;  ++ii) {SubdiagClear(ii+block);  subdiags[ii+block][0] = subdiags[ii][0];}
  dim += nn;

  return;
}



// BuildVecs
// Given a set of 'numvecs' x 'lenvecs' 'coefficients' and a place ('toBuild') to put the vectors it builds, this builds the vectors.
// It assumes either that the Lanczos procedure has not yet been started or that all necessary vectors have been stored (scratch==NULL).
// (it would be possible to generate those vectors not yet stored, but this might be a useless generalization.)
// Since it is now possible to add multiple vectors to the Lanczos space in one iteration, the 'modulus' option 
// (which must be <= block) says how many to generate at once.
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::BuildVecs(LancVec **toBuild, LancNum **coefficients, int numvecs, int lenvecs, int modulus)
	{
	if (dim == 0)
		{int ii = 0;
		 for (int jj=0; jj<numvecs; jj++) {                toBuild[jj]->LinCombo(0,*start[0],   coefficients[jj][ii],*start[ii]);}
		 for (++ii; ii<block && ii<lenvecs; ii++)
			{        for (int jj=0; jj<numvecs; jj++) {toBuild[jj]->LinCombo(1,*toBuild[jj],coefficients[jj][ii],*start[ii]);}}
		 LancVec** current_set = vectors + block;
		 while (ii < lenvecs)
			{int numnew = std::min(modulus,lenvecs-ii);
			 AddVectors(numnew);
			 for (int kk=block-numnew; kk<block; kk++, ii++)
			 	{for (int jj=0; jj<numvecs; jj++) {toBuild[jj]->LinCombo(1,*toBuild[jj],coefficients[jj][ii],*current_set[kk]);}}}}
	else
		{int ii = 0;
		 for (int jj=0; jj<numvecs; jj++) {                toBuild[jj]->LinCombo(0,*start[0],   coefficients[jj][ii],*start[ii]);}
		 for (++ii; ii<lenvecs; ii++)
			{        for (int jj=0; jj<numvecs; jj++) {toBuild[jj]->LinCombo(1,*toBuild[jj],coefficients[jj][ii],*start[ii]);}}}
	return;
	}





// // // // // // // // // // // // // // // // Private Funtions // // // // // // // // // // // // // // // // // //



// Initialize
// This populates an array whose elements point to the vectors we are currently working with, ie, those vectors
// which will have the matrix acting on them (the middle third), the locations to store the result of this
// (the last third) and those vectors whose only purpose is to be Gram-Schmidt-ed against (the first third).
// Mostly we will proceed by simply shifting the pointer cyclically after any call to AddVectors, unless the
// Lanczos basis is being saved, in which case "fresh" vectors are used to store the results of these operations.
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::Initialize(bool lanczosbasis)
	{
	dim = 0;								// The number of new vectors generated so far
	
	if (vectors) {delete [] vectors;}  vectors = new LancVec* [3*block];	// Set up array for managing memory shuffling
	int ii = 0;
	for (int jj=0; jj<block; ii++, jj++) {vectors[ii] =        NULL;}
	for (int jj=0; jj<block; ii++, jj++) {vectors[ii] =   start[jj];}
	for (int jj=0; jj<block; ii++, jj++) {vectors[ii] = scratch[jj];}
	for (int jj=0; jj<block;       jj++) {         SubdiagClear(jj);}	// Our 'block'x'block' approx to the projection so far

	if (lanczosbasis) {scratch = NULL;}					// This pointer value is used as a 'lanczosbasis' signal
	return;
	}



// SubdiagClear
// Private utility to set a column of the projection to zero
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::SubdiagClear(int ii) {for (int jj=0; jj<=block; jj++) {subdiags[ii][jj]=0;}  return;}



// ShiftPointers
// Private utility to shift the pointers after nn new vectors are made
//
template <typename LancNum, class LancVec, class LancMat>
void LanczosEngine<LancNum,LancVec,LancMat>::ShiftPointers(int nn)
	{

	for (int ii=0; ii<3*block-nn; ii++) {vectors[ii] = vectors[ii+nn];}

	if (scratch) {for (int ii=3*block-nn; ii<3*block; ii++) {vectors[ii] = vectors[ii-2*block];}}
	else         {for (int ii=3*block-nn; ii<3*block; ii++) {vectors[ii] = start[dim+2*block+ii];}}		// assumes dim not yet incremented

	return;
	}



}	// End namespace LanczosProjection

#endif	// End #ifdef _LANCZOSENGINEH_

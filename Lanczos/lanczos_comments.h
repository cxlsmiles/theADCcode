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

/*****     A (band-)Lanczos-projection-based-algorithm Suite for Eigenvalue/vector Refinement or Propagation Methods     ***************
****************************************************************************************************************************************

The 'Lanczos' class defined here contains a number of utilities associated with the output of a private workhorse class for Lanczos 
projection ('LanczosEngine').  The idea is that the desired Lanczos-based algorithm may be built transparently on top of these utilities,
of which, more can always be added.  In order to leave the underlying engine most flexible and uncluttered (also faster and predicatable), 
its calling protocol is simple (but somewhat tedious to work with), and unprotected (in that the engine does less error checking).  It is 
the job of these 'Lanczos' utilities to know how to use the 'LanczosEngine' class, and use it safely, performing error checking first, and 
doing memory management outside of the engine.  This "comment" serves as the introductory documentation.

Contents:
1. What is Lanczos projection, abstractly?
2. How to use the 'Lanczos' class
   -subpoint:  What properties the template data types must have
3. How to interact with the 'LanczosEngine' class and its output (for building more 'Lanczos' utilities)


	1. What is Lanczos projection, abstractly?

The common task among all algorithms is to build a projection of a Hermitian linear operator (henceforth refered to as a matrix, since 
computers are limited to discrete algebras) which is either too big to be built entirely, or too big to diagonalize or otherwise manipulate
as desired.  This section describes the properties of such a projection on a level which is more abstract than code protocol, but addresses
some of the properties of the projection that the user will have control over.

A matrix supplied by the user is projected into the Krylov subspace defined by that matrix and a starting vector (also supplied by the 
user), or the direct sum of such Krylov subspaces, if a set (often refered to as a block) of initial vectors is given.  If the block size 
is one, a simple Lanczos projection is performed (as opposed to block/band Lanczos, as described below).  For the uninitiated, the Krylov 
subspace of order N, defined by a vector and a matrix, is the span of the set of N+1 vectors generated by acting the matrix on the vector 
0 through N times.  The internal representation of the vectors is left to the user, as long as they form a linear inner-product space, with
a positive-definite metric.  Of course certain public functions for taking linear combinations and inner products must be defined, and this
is discussed later.  Likewise, the representation of the matrix is left to the user, and it need not be a two-dimensional array of matrix
elements.  In fact, only its action on a vector is ever needed.  This is particularly valuable when the matrix would be too large to store.
If matrix is not Hermitian the algorithms here need deep modification, since the projected matrix is upper/lower triangular (plus some 
bands), rather than band diagonal, as will be discussed.

A starting vector, which defines a Krylov subspace, must be normalized, and if a block of such vectors are used, they need to form an
orthonormal set.  The (direct sum of) Krylov space(s) is then spanned by the orthonormal Lanczos vectors, which are generated 
systematically by this code and the previous Lanczos vectors.  The representation of the projection of the matrix built here is indeed 
that of a two-dimensional array of matrix elements, in the basis of these Lanczos vectors.  The starting vectors (always refered to in the
plural from now on, for convenience) are the first members of this Lanczos basis by definition.  The further Lanczos vectors beyond the 
starting block are defined as in the article by Weikert, Meyer and Cederbaum [J. Chem. Phys., v.104 p.7122 (1996)] (i.e. the manner in 
which the subsequent Krylov vectors are orthonormalized to make Lanczos vectors).  

This code operates a band Lanczos algorithm, which is the same as a block Lanczos algorithm for projections whose dimensionality are an 
integer multiple of the block size, except for the exact definition of the orthonormal Lanczos vectors spanning the Krylov space 
(implicitly meaning the direct sum of Krylov spaces, from now on, when relevant).  Unlike the block algorithm, the dimensionality of the 
projection need not be an integer multiple of the block size, but if it is not (dimensionality is specified by the user in a variety of 
manners), then the Krylov subspaces built from the starting vectors with lower index (as given in input) will be of one order higher than
the others.  When block size is one, the block and band algorithms are identical in every respect.  In a band algorithm, the projection of
a Hermitian matrix is band diagonal with a band width such that it has B nonzero superdiagonal and B nonzero subdiagonal elements in each
column, where B is the number of vectors in the starting block.  This is a more compact and uniform representation than a block algorithm
which results in a block-tridiagonal matrix of BxB blocks (simple tridiagonal for simple single-vector Lanczos in each case).  The 
projection generated is Hermitian, so there are only B+1 elements that need to be saved from each column of the projection in the 
band-diagonal matrix (the diagonal and B offdiagonals).  

One important convention needs to be made clear here.  In short, for any DxD projection, the upper-left (D-B)x(D-B) submatrix is exact, 
along with all other elements (some off-diagonals, discussed below) that are not in the lower-right BxB submatrix, where B is the size of 
the starting block.  This final block has approximate diagonals and zero as off-diagonals.  The common approximation of setting the final
B diagonals equal to the preceding B diagonals is used.  By convention, if D=B, then we just have a BxB matrix of zeros as the approximate
projection, with no exact part yet.  (Essentially, this is the initialization point where the starting block has been defined, but no work
has been done.)  At every iteration, these approximate and zeroed values are set as such in the projection, but one can choose not to use
them by simply ignoring the final B columns of the output projection (and some trailing members of the final B columns of what is left
after that) and taking only the exact part (of smaller dimension).  It is worth noting that the exact off-diagonals falling outside of the
exact and approximated submatrices are useful for computing eigenvector residuals, exactly (from diagonalization of the exact part), along
with some other convergence metrics for other methods.  To understand why this is, one needs to know that, due to the details of a Lanczos
algorithm, for any exact representation of the matrix in a basis of Lanczos vectors up to a Krylov given order, a number of further 
offdiagonal elements are computed outside of this exact projection.  These correspond to elements between vectors inside the space where 
the projection is exact and vectors of higher order.  However, these elements come without knowing the corresponding diagonal elements 
(or some elements less removed from the diagonal) between the higher order vectors.  Attempting to compute these diagonal (or 
less-offdiagonal) elements, to increase the size of the exact projection, results in "accidental" computation of further such offdiagonals
outside of the "new" exact part (for which the diagonals, etc., are again not known).  Therefore, this code employs the common 
approximation for the unknown values automatically, so as to not discard the information that we do know.  


	2. How to use the 'Lanczos' class

The class template depends on three data types (classes), a scalar class, a vector class and a matrix class, whose template parameter 
names are LancNum, LancVec and LancMat, respectively.  The expected properties of these classes is the last part of this section.  What 
follows until then is a rough guide to the way the Lanczos class operates.  Details on the exact arguments of the constructor and member 
functions are to be found near the relevant parts of the code.  The code centers around an object of type LancMat, which is some 
representation of the matrix to be projected, and this object is the only mandatory argument of the constructor.  Exactly how this matrix
is accessed (its action on a vector) by the code here is addressed in the class properties information, to follow at the end of this 
section.

The main job of the constructor is simply to set up a bunch of internal bookkeeping variables.  One aspect of the internal design of the 
Lanczos object needs to be addressed, because it affects how the user interacts with it, and we will need to refer to it later.  In the 
design of this class, one encounters the common problem that the amount of information stored depends on how many iterations of some task
are performed, in this case, the building of further projection columns.  One then needs a mechanism to keep track of these allocation
pointers in an orderly fashion.  Linked lists have their own difficulties in terms of memory overhead per "array element" and overhead in
terms of lines of code (and CPU time) to access elements.  This code therefore opts for simply allocating a couple of very long private
arrays of pointers in the constructor, the elements of which are to be assigned as needed.  One needs one element of one array for each 
column of the projection to be built (and potentially one element of the other for each Lanczos vector, if they are saved).  This is a 
relatively safe strategy, since Lanczos projections are meant to be small (to be useful, and for stability), and so static arrays of a 
size much bigger than necessary contribute only a small memory overhead; the user may override the default size of these arrays.  Many 
Lanczos algorithms can also be restarted using information from a run which has exhausted its pointers.  

A second important thing to note is that this algorithm is built to be stingy with memory.  This mostly refers to not storing any more 
vectors than necessary, since it is assumed that all other information is small by comparison.  The smaller storage requirements of 
Lanczos-based algorithms is one of their strongest advantages.  This is the reason why some input is overwritten, as explained below.


The general flow is as follows:

a) Declare or allocate an object of the 'Lanczos' class, specifying what data types/classes will be used for scalars, vectors and matrices,
   and specify the matrix to be projected.
b) Allocate and build the orthonormal starting vectors outside of the Lanczos class, and pass in the pointers using the member functions
   AddStartVec() or AddStartBlock().  For what is believed to be reasons of good code design, the responsibility for freeing these vectors
   lies with the calling function that allocated them.  However, the user should be aware that these vectors will be overwritten.  The
   Lanczos routines here use both these vectors and some internal space (for which it is responsible) as scratch space.  The exception to
   this rule is when the Lanczos basis is saved, in which case the starting vectors are saved as well, since they are members of this 
   basis.  The calling function should retain a copy of the start-vector pointers, both to free the vectors when the time comes
   and also because it may be necessary to remake the starting vectors inside of these allocations, in cases where Lanczos will be run 
   twice (building some desired linear combination(s) of the Lanczos vectors in the second run, based on some result of the first).  These
   functions to add start vectors may be called as many times as desired to add more vectors to the start block, up until Iterate() is 
   called the first time.  They automatically update the variable that holds the size of the starting block.
c) Call Iterate() to allocate additional internal space for the columns of the projection, and generate said information.  This function
   adds an entire block of the matrix, the dimension of the starting block.  These columns are pointed to by the successive elements of
   one of the internal arrays above.  If storing of the Lanczos basis is desired, then additional vectors are allocated, pointed to by the
   members of the other internal array.  Otherwise only the first block of elements of this array are used to hold the start vectors 
   (which are then overwritten).  Iterate() can basically be called as many times as the user desires, until the projection has somehow
   converged (or the number of internal pointers is exhausted).  After the first call to Iterate(), an approximate BxB matrix is built,
   where B is the number of vectors in the starting block.  The matrix is all zeros, as explained in the section 1.  After the second 
   call, an approximate 2Bx2B matrix is built, where the upper-left BxB block is exact, along with the upper-right and lower-left BxB
   off-diagonal blocks, and the lower-right BxB block is approximated.  Subsequent calls follow the pattern outlined in the section 1,
   building successively larger exact blocks, always with a trailing BxB approximate block.
d) The next steps depend very much on the user.  For example, there is a Diagonalize() routine which diagonalizes the most recent 
   projection, and a DiagConverged() routine that checks the residuals of the eigenvectors which have the largest overlaps with the start
   block (assumes only the exact part of the matrix was diagonalized, which the user can control).  Also the function GetSubdiags() will
   return the pointer to the internal array that holds the projection.  This output is stored as D (B+1)-dimensional arrays, where D is
   the dimensionality of the projection and B is the number of vectors in the start block.  The 0th element of each array is a diagonal 
   element and the subsequent members are the B nonzero subdiagonals for that column, which define the B nonzero superdiagonals by 
   Hermitian symmetry.  Taking the complex conjugate of each element converts it to storage by rows, instead of columns.
e) A very common operation will be to build those eigenvectors of the projection with the largest overlap with the starting vectors, in
   the representation of the original matrix.  The function BuildEigenvecs() does exactly this.  One allocates the vectors to be filled
   with these eigenvectors outside of the Lanczos class, and passes in the pointers.  The vectors are filled and a pointer to an internal
   array with the respective eigenvalues (both arranged in ascending order of the eigenvalues) is returned.  Since this routine will
   have to run Lanczos again (if the Lanczos vectors were not saved), the guess vectors will need to be built again (in place) before
   this is called.  Usually the Lanczos vectors are too expensive to save, but if they were, the guess vectors were not overwritten and
   Lanczos does not need to be run again, and BuildVectors() will take appropriate action.  No new memory is allocated internally, since
   everything that would be needed for a second run has already been allocated during the first one.
f) The function CleanUp() is exactly what is called by the destructor.  It may be called prematurely to release memory allocated inside
   this class, and everything will be destroyed that was not allocated outside of the class (and none of that will be released).

The data types/classes are expected to have the following properties:

scalar class: template parameter = 'LancNum'
The flexibility in the scalar class allows for computation in single or double (or arbitrary) precision with real or complex data.  Most 
generally, the operations which are defined for single or double precision floating-point types should be valid for this type, but the 
following are minimal requirements for the scalar class (if type float or double is used, these are defined by default):
	unary  -       : returns the negative of the operand as the same data type.
	unary  +       : returns the complex conjugate of the operand as the same data type.
	unary (double) : This should return a double which is the real part of the operand.
	binary =       : in addition to the default copy assignment, assignment from any real native type (at least double and integer) 
	                 should be defined as assigning the real part (imaginary part is zero) and this should return the value of the 
	                 post-operation left-hand operand as the same data type.
	binary += , -= : This should increment the left-hand operand by the value of the right-hand operand (or its negative, 
	                 respectively), which is of the same data type, and return the post-operation left-hand operand.
	binary *       : This should return the product of the left and right operands, both of the same data type as the return value.

vector class: template parameter = 'LancVec'
The vector class used needs to have the following operations and functions defined:
	copy constructor : It is assumed that a vector object may exist even if no storage (core or disk) is associated with it.  However,
	              it is not assumed that such an object is compact, even if no storage is associated with it (yet).  Therefore, the 
	              underlying algorithms operate largely in terms of arrays of pointers to vectors, so that additional vector objects 
	              may be allocated only as needed.  For such an allocation, a copy constructor is necessary, so that the attributes of
	              one of the starting vectors may simply be copied.  However, the new vector must not point to the storage of the 
	              copied vector.  Upon copy construction, only such parameters (like dimensionality) that indicate that the vectors 
	              belong to the same space should be copied, but the storage associated with it needs to be unique.  The new vector 
	              need not represent a copy of the copied vector in the space, i.e. the amplitudes to not need to be copied (It would 
	              use extra time, but not hurt anything.).  Therefore the copy constructor can either allocate new storage (perhaps 
	              leaving it empty), or create a vector that points to no storage, if storage is allocated by the vector class, as 
	              needed during manipulations.
	binary   || : returns the l-norm of the vector as a real type (double, float, etc) from which type LancNum can take assignment.
	              The second operand is an integer which specifies l.  Only l=2 needs to be implemented, and it should return what is
	              usually called the norm-squared of the vector.
	binary   *  : returns the inner product of two vectors as any type from which type LancNum can take assignment (probably LancNum).
	              The inner product is defined as the matrix (dot) product of the Hermitian conjugate of one vector with another.
	binary   /= : scales the vector (first operand) by the inverse of the second operand, which is of any type from which LancNum
	              can take an assignment (at least LancNum itself and double).
	function LinCombo(number,const LancVec,number,const LancVec) : assigns the vector to be the linear combination of the two vectors 
	              that are given as arguments, each weighted by the scalar argument that immediately precedes it.  The implementation 
	              should allow that the vector to be assigned to is also one of the arguments (specifically, the first of the two 
	              vectors).  [It is recommended that this function be implemented as a template, so that the types 'number' maybe any 
	              scalar type, including integers.  (Specifically, this function will be called with the first coefficient of type int
	              and the second of type LancNum.)  This can easily be done by wrapping a function where the scalar types are known to 
	              be of type LancNum, and and passing the wrapped function LancNum parameters which have been assigned from these 
	              parameters of LinCombo(), since LancNum can be assigned from any scalar type (as specified above).  For efficiency,
	              the implementation may want to check whether the first two arguments have the value '1' and '*this', respectively.]

matrix class: template parameter = 'LancMat'
The matrix class used needs to have the following operator defined:
	unary (LancVec, const LancVec) : This operator only reads the second argument and then writes the product of the matrix onto that
	              vector into the first argument.  [This may be implemented as (LancVec&, const LancVec&), especially if there is no 
	              LancVec default constructor.]


	3. How to interact the 'LanczosEngine' class and its output

At the time of construction, the 'LanczosEngine' object expects all of the information necessary to do its work to be in place.  The first
argument, of type LancMat*, is the pointer to the matrix to be projected.  The second argument gives the number of vectors in the starting
block, which is a constant for the engine.  An array of LancVec* pointers to the starting vectors is the third argument, and the forth 
argument is an array of pointers to LancNum arrays in which to put the subdiagonals of the projected matrix.

One of the possible outputs of the 'LanczosEngine' object is the Lanczos vector basis.  This can potentially cost a lot of storage space 
and is not always necessary.  Whether to store this basis is a decision of a given algorithm.  This is indicated by the final Boolean
argument.  Since the starting vectors are the first members of the Lanzos basis, it is assumed that pointers to the further vectors will 
be stored in the array of start vectors, and in this case, that array will need to be longer than the start block size.  As it the engine
is presently used, it is always the same length as the array for the projection, but the later elements will not be referenced if the 
Lanczos basis is not saved.  If the Lanczos basis is not being stored, it is assumed that storage space is to be conserved, and the 
starting vectors are overwritten eventually for scratch space (the unused pointers are a comparitive minimal overhead).

Immediately after the 'LanczosEngine' constructor has been called, and before this engine has been asked to generate any further Lanczos 
vectors via the functions AddVector() or AddBlock() (no matrix-vector operations have been done yet) the first B columns of the projection 
are already filled in, where B is the number of vectors in the starting block, though they are just zero, since we know absolutely nothing
at this point.  This means that the first B entries of the array for the projection need to point to arrays of B+1 LancNum elements 
that have already been allocated.  Likewise, the first B elements of the array of start vectors need to point to the start vectors
(as written above).  However, the later elements of both of these arrays have not yet been referenced, and one may withhold assigning these
to point to allocations of further projection columns and vectors until one is prepared to add that column, and the latter allocation
is unnecessary, unless one is saving the Lanczos basis.

The primary output is are the B+1 nonzero subdiagonal elements of the matrix projection, which specify the B+1 nonzero superdiagonals by
Hermitian symmetry.  They are in the same format as described in section 2, step d).  Indeed, it is the same array.

****************************************************************************************************************************************
****************************************************************************************************************************************/

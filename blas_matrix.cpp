#include "blas_matrix.hpp"
#include "integral_blocks.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


namespace lapack {
  extern "C" {

#if defined  INTEL
#include <mkl_blas.h>
#include <mkl_lapack.h>
#elif defined PGI
#include <acml.h>
#elif defined GCC
    void dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
    double dnrm2_(int *n, double *x, int *incx);
    double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
    void dscal_(int *n, double *alpha, double *x, int *incx);
    void daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy);
    void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
                int *lda, double *x, int *incx,	double *beta, double *y, int *incy);
    void  dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
                 double *a, int *lda, double *b, int *ldb, double *beta, double *c, int	*ldc);
    void dsymv_(char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx,
                double *beta, double *y, int *incy);
    void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, 
                double *work, int *lwork, int *info);
    void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
    int ilaenv_(int *ispec, char *name__, char *opts, int *n1, int *n2, int *n3, int *n4);
    void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
#endif

  }
}

// Turns on some safety checks (for debugging)
#define SAFE
#undef  SAFE

//Prints the matrix contents
void Blas_matrix::print() 
{
  
#define COLUMNS_PER_LINE 8
#define PRINT_ALL
#define PRECISION 8

  // Check if the number of columns is not smaller than 
  // the columns-per-line constant.
  int col_per_line = (COLUMNS_PER_LINE > cols()) 
    ? cols() : COLUMNS_PER_LINE;
  
  int current_col = 0;
  int column_group = 0;
  // Print while the final column is reached
  while (current_col < cols()) {
    // Determine the starting and ending column 
    // for the current column group that is to be printed.
    int start_col = column_group++ * col_per_line;
    int end_col = column_group * col_per_line - 1;
    // Check if the end has been reached
    end_col = (end_col >= cols()) ? cols() - 1 : end_col;
    
    // Print the column numbers
    cout << "\n Columns " << start_col << " through " << end_col << endl;
    for(current_col = start_col; current_col <= end_col; current_col++)
      cout << '\t' << setw(PRECISION+2) << current_col;
    cout << endl;
    
    // Print all rows
    for(int current_row = 0; current_row < rows(); current_row++) {
      
      bool printed = false; 
      
      // Print the data
      for(current_col = start_col; current_col <= end_col; current_col++) {
#ifndef PRINT_ALL
	//Print only the lower triangle
	if (current_col > current_row) break;
#endif
	if (!printed)
	  cout << setw(3) << current_row;

	printed = true;
	

	cout << '\t' << setiosflags(ios::fixed) << setprecision(PRECISION) 
	     << (*rec_)(current_row, current_col);
      }

      if (printed)
	cout << endl;
      
    }
  }
}



// Copies some external memory if given any
Blas_matrix::Blas_matrix(unsigned int rows, unsigned int cols, double *mat)
{
  rec_ = new Rectangular_matrix<double>(rows, cols);
  
  if (mat) {
    int inc = 1, length = size();
#if defined  INTEL
    lapack::dcopy(&length, mat, &inc, &(*rec_)(0,0), &inc);
#elif defined PGI
     lapack::dcopy(length, mat, inc, &(*rec_)(0,0), inc);
#elif defined  GCC
    lapack::dcopy_(&length, mat, &inc, &(*rec_)(0,0), &inc);
#endif
  }
  //else 
  //  dcopy(0.);
 
}

// A requirement by the Lanczos template: a copy constructor
Blas_matrix::Blas_matrix(const Blas_matrix &bx)
{
  rec_ = new Rectangular_matrix<double>(bx.rows(), bx.cols());
  
  int inc = 1, length = size();
#if defined  INTEL
  lapack::dcopy( &length, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::dcopy( length, &bx(0,0), inc, &(*rec_)(0,0), inc);
#elif defined  GCC
  lapack::dcopy_( &length, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#endif 
  
}

// A dcopy wrapper
void Blas_matrix::dcopy(Blas_matrix &bx)
{

#ifdef SAFE
  if (size() < bx.size()) 
    throw string("Blas_matrix::dcopy:source larger than the copy");
#endif
  
  int inc = 1, size = bx.size();
#if defined  INTEL
  lapack::dcopy( &size, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::dcopy( size, &bx(0,0), inc, &(*rec_)(0,0), inc);
#elif defined  GCC
  lapack::dcopy_( &size, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#endif
}

// A dcopy wrapper
void Blas_matrix::dcopy(Submatrix &bx)
{

#ifdef SAFE
  if((bx.row_slice != rows()) || ( bx.col_slice != cols()))
    throw std::string("Blasmatrix::dcopy=(): Different matrix dimensions.");
#endif
  int inc = 1;
  int size = bx.row_slice;
  for (unsigned col = 0; col < bx.col_slice; col++)
#if defined  INTEL
    lapack::dcopy(&size, &bx.matrix(bx.row_start, col + bx.col_start), &inc, 
    	&(*rec_)(0,col), &inc); 
#elif defined PGI 
    lapack::dcopy(size, &bx.matrix(bx.row_start, col + bx.col_start), inc, 
		  &(*rec_)(0,col), inc); 
#elif defined  GCC
    lapack::dcopy_(&size, &bx.matrix(bx.row_start, col + bx.col_start), &inc, 
    	&(*rec_)(0,col), &inc); 
#endif
}



// A dcopy wrapper
void Blas_matrix::dcopy(double alpha)
{
  int inc = 1, noinc = 0, length = size();
#if defined  INTEL
  lapack::dcopy(&length, &alpha, &noinc, &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::dcopy(length, &alpha, noinc, &(*rec_)(0,0), inc);
#elif defined  GCC
  lapack::dcopy_(&length, &alpha, &noinc, &(*rec_)(0,0), &inc);
#endif
  
}

//A dcopy wrapper overloading the assignment operator
Blas_matrix& Blas_matrix::operator=(Blas_matrix &bx)
{
  this->dcopy(bx);
  return *this;
}


//A dcopy wrapper: copies a real number to all elements
Blas_matrix& Blas_matrix::operator=(double alpha)
{
  this->dcopy(alpha);
  return *this;
}



//Required by the Lanczos template: a dnrm2 wrapper.
// Note that the result is squared!!!
double Blas_matrix::operator||(unsigned int l) 
{
#ifdef SAFE
  if (l != 2)
    throw std::string("Blas_matrix::operator||: Only Euclidean norm implemented.");
#endif
  
  int inc = 1, length = size();
#if defined  INTEL
  double norm = lapack::dnrm2(&length, &(*rec_)(0,0), &inc);
#elif defined PGI
  double norm = lapack::dnrm2(length, &(*rec_)(0,0), inc);
#elif defined  GCC
  double norm = lapack::dnrm2_(&length, &(*rec_)(0,0), &inc);
#endif
  return norm * norm;
}
  

// A dscal wrapper required by Lanczos
 Blas_matrix& Blas_matrix::operator/=(double alpha)
{
#ifdef SAFE
  if (!alpha) 
    throw std::string("Blas_matrix::operator/=: Division by zero.");
#endif
  
  alpha = 1. / alpha;
  
  this->dscal(alpha);
  return *this;
}



//A ddot wrapper
 double Blas_matrix::ddot(Blas_matrix &bx) 
{
#ifdef SAFE
  if (size() != bx.size())
    throw std::string("Blas_matrix::ddot:  Different vector sizes.");
#endif

  int inc = 1, length = size();
#if defined  INTEL
  return lapack::ddot(&length, &(*rec_)(0,0), &inc, &bx(0,0), &inc);
#elif defined PGI
  return lapack::ddot(length, &(*rec_)(0,0), inc, &bx(0,0), inc);
#elif defined  GCC
  return lapack::ddot_(&length, &(*rec_)(0,0), &inc, &bx(0,0), &inc);
#endif  
}


// A ddot wrapper
double Blas_matrix::ddot(const Submatrix &x) 
{
#ifdef SAFE
  if((x.col_slice != 1) || (cols() != 1))
    throw std::string("Blas_matrix::ddot(): Term not a vector.");
  if(x.row_slice != size())
    throw std::string("Blas_matrix::ddot():  Different vector sizes.");
#endif
  
  int inc = 1;
  int size = x.row_slice;
#if defined  INTEL
   return lapack::ddot(&size, &(*rec_)(0,0), &inc, 
 		    &x.matrix(x.row_start, x.col_start), &inc);
#elif defined PGI
  return lapack::ddot(size, &(*rec_)(0,0), inc, 
		    &x.matrix(x.row_start, x.col_start), inc);
#elif defined  GCC
   return lapack::ddot_(&size, &(*rec_)(0,0), &inc, 
 		    &x.matrix(x.row_start, x.col_start), &inc);
#endif  
}


//A ddot wrapper overloading the multiplication operator
double Blas_matrix::operator*(Blas_matrix &bx) 
{
  return this->ddot(bx);  
}

//A ddot wrapper overloading the multiplication operator
double Blas_matrix::operator*(const Submatrix &bx) 
{
  return this->ddot(bx);  
}



// A dscal wrapper
 void Blas_matrix::dscal(double alpha)
{
  int inc = 1, length = size();
#if defined  INTEL
  lapack::dscal(&length, &alpha, &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::dscal(length, alpha, &(*rec_)(0,0), inc);
#elif defined  GCC
  lapack::dscal_(&length, &alpha, &(*rec_)(0,0), &inc);
#endif
}

// A dscal wrapper overloading the assignment operator
Blas_matrix& Blas_matrix::operator*=(double alpha)
{
  this->dscal(alpha);
  return *this;
}

// A daxpy wrapper
void Blas_matrix::daxpy(double alpha, Blas_matrix &bx)
{

#ifdef SAFE
  if (size() != bx.size())
    throw std::string("Blas_matrix::daxpy:  Different vector sizes.");
#endif  

  int inc = 1, length = size();
#if defined  INTEL
  lapack::daxpy(&length, &alpha, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::daxpy(length, alpha, &bx(0,0), inc, &(*rec_)(0,0), inc);
#elif defined  GCC
  lapack::daxpy_(&length, &alpha, &bx(0,0), &inc, &(*rec_)(0,0), &inc);
#endif  
}




// This set of daxpy methods gets as an argument a
// specific subarea of the integral tables whose memory structure is different in each case.
// The daxpy method is overloaded for each particular type of "argument".

// Here a polymorphic call is used to pass a reference to the current Blas_matrix.
// The caller will then use the passed reference to call one of the overloaded
// daxpy versions bellow. (see mat_env.hpp)
void Blas_matrix::daxpy(double alpha, Daxpy_argument* b)
{
  b->get_argument(alpha, *this);
}



// Adds a rectangular matrix to the local matrix.
// This rectangular matrix is expected to be contiguous 
// just along its columns. Thus daxpy is called for
// each column.
void Blas_matrix::daxpy(double alpha, Rectangular_argument &r)
{
#ifdef SAFE
  if ((rows() != r.rows()) || (cols() != r.cols()))
    throw std::string("Blas_matrix::daxpy: rec Different vector sizes.");
#endif  

  int inc = 1, length = rows();
  
  for(unsigned int col = 0; col < cols(); col++)
#if defined  INTEL
    lapack::daxpy(&length, &alpha, &r(0,col), &inc, &(*rec_)(0,col), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &r(0,col), inc, &(*rec_)(0,col), inc);
#elif defined  GCC
    lapack::daxpy_(&length, &alpha, &r(0,col), &inc, &(*rec_)(0,col), &inc);
#endif
}



// Adds a transposed matrix to the local matrix.
// The transposed matrix is contiguous along its rows.
// daxpy is called for each row and the elements are
// added to the rows of the local matrix by using a
// step equal to the column size
void Blas_matrix::daxpy(double alpha, Transposed_argument &t)
{
#ifdef SAFE
  if ((rows() != t.rows()) || (cols() != t.cols()))
    throw std::string("Blas_matrix::daxpy: trans Different vector sizes.");
#endif  

  int inc = 1, inc_row = rows(), length = cols();
  
  for(unsigned int row = 0; row < rows(); row++)
#if defined  INTEL
    lapack::daxpy(&length, &alpha, &t(row,0), &inc, &(*rec_)(row,0), &inc_row);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(row,0), inc, &(*rec_)(row,0), inc_row);
#elif defined  GCC
    lapack::daxpy_(&length, &alpha, &t(row,0), &inc, &(*rec_)(row,0), &inc_row);
#endif 
}


// The argument here is a triangular matrix whose columns are contiguous chunks.
// Each of its columns is first added to the corresponding column of the local matrix
// (from the beginning of the column up to the relevant diagonal element)
// and then to the row of the local matrix with the same number 
// (up to the element preceding the diagonal element)
void Blas_matrix::daxpy(double alpha, Upper_triangular_argument& t)
{
#ifdef SAFE
  if ((rows() != t.rows()) || (cols() != t.cols()))
    throw std::string("Blas_matrix::daxpy:  Different vector sizes.");
#endif  

  int inc = 1, inc_row = rows();

  for(unsigned int col = 0; col < cols(); col++) {
    int length = col + 1;
#if defined  INTEL    
    lapack::daxpy(&length, &alpha, &t(0,col), &inc, &(*rec_)(0,col), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(0,col), inc, &(*rec_)(0,col), inc);
#elif defined  GCC    
    lapack::daxpy_(&length, &alpha, &t(0,col), &inc, &(*rec_)(0,col), &inc);
#endif
    length--;
#if defined  INTEL
    lapack::daxpy(&length, &alpha, &t(0,col), &inc, &(*rec_)(col,0), &inc_row);  
#elif defined PGI    
   lapack::daxpy(length, alpha, &t(0,col), inc, &(*rec_)(col,0), inc_row);    
#elif defined  GCC
    lapack::daxpy_(&length, &alpha, &t(0,col), &inc, &(*rec_)(col,0), &inc_row);  
#endif  
  }  

}


// Similar to above, however here the contiguous columns are added to
// lower triangale part of the columns, and to the rows in the upper triangle
// to the corresponding rows/columns of the local matrix
void Blas_matrix::daxpy(double alpha, Lower_triangular_argument& t)
{

#ifdef SAFE
  if ((rows() != t.rows()) || (cols() != t.cols()))
    throw std::string("Blas_matrix::daxpy:  Different vector sizes.");
#endif  
  
  int inc = 1, inc_row = rows();

  for(unsigned int row = 0; row < rows(); row++) {
    
    int length = rows() - row;
#if defined  INTEL    
    lapack::daxpy(&length, &alpha, &t(row,row), &inc, &(*rec_)(row,row), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(row,row), inc, &(*rec_)(row,row), inc);
#elif defined  GCC   
    lapack::daxpy_(&length, &alpha, &t(row,row), &inc, &(*rec_)(row,row), &inc);
#endif
    length--;
    
    if (length)
#if defined  INTEL
    lapack::daxpy(&length, &alpha, &t(row+1,row), &inc, &(*rec_)(row,row+1), &inc_row);    
#elif defined PGI
    lapack::daxpy(length, alpha, &t(row+1,row), inc, &(*rec_)(row,row+1), inc_row);    
#elif defined  GCC
    lapack::daxpy_(&length, &alpha, &t(row+1,row), &inc, &(*rec_)(row,row+1), &inc_row);    
#endif

  }
  
}


// Adds a vector to the diagonal of the matrix
void Blas_matrix::add_diag(Blas_matrix &bv)
{
#ifdef SAFE
  if (rows() != cols())
    throw std::string("Blas_matrix::add_diag: Not a square matrix.");
  if (rows() != bv.size())
    throw std::string("Blas_matrix::add_diag: Vector length different from the diagonal length.");
#endif

  double alpha = 1.;
  int inc = 1, diag_inc = rows() + 1, size = rows();

#if defined  INTEL 
  lapack::daxpy(&size, &alpha, &bv(0,0), &inc, &(*rec_)(0,0), &diag_inc);
#elif defined PGI
  lapack::daxpy(size, alpha, &bv(0,0), inc, &(*rec_)(0,0), diag_inc);
#elif defined  GCC 
  lapack::daxpy_(&size, &alpha, &bv(0,0), &inc, &(*rec_)(0,0), &diag_inc);
#endif
  
}


// Adds a number to the diagonal of the matrix
void Blas_matrix::add_diag(double alpha)
{
#ifdef SAFE
  if (rows() != cols())
    throw std::string("Blas_matrix::add_diag: Not a square matrix.");
#endif

  double one = 1.;  
  int no_inc = 0, diag_inc = rows() + 1, size = rows();

#if defined  INTEL
  lapack::daxpy(&size, &one, &alpha, &no_inc, &(*rec_)(0,0), &diag_inc);
#elif defined PGI
  lapack::daxpy(size, one, &alpha, no_inc, &(*rec_)(0,0), diag_inc);
#elif defined  GCC
  lapack::daxpy_(&size, &one, &alpha, &no_inc, &(*rec_)(0,0), &diag_inc);
#endif  
}



// A dgemv wrapper: Performs y = A.x + y, where x and y are vectors
void Blas_matrix::dgemv(char transa, Blas_matrix &A, Blas_matrix &x)
{
#ifdef SAFE
  if((cols() != 1) || (x.cols() != 1))
    throw std::string("Blas_matrix::dgemv(): Term not a vector.");
  
  if(
     ((transa == 'N') && ((rows() != A.rows()) || 
			  (A.cols() != x.rows())))
     ||
     ((transa == 'T') && ((rows() != A.cols()) || 
			  (A.rows() != x.rows())))
     )
    throw std::string("Blas_matrix::dgemv(): Different matrix sizes.");
#endif  
  
  double one = 1.;  
  int rows = A.rows(), cols = A.cols(), inc = 1, lda = rows;
#if defined  INTEL
   lapack::dgemv(&transa, &rows, &cols, &one, &A(0,0), &lda, 
 		 &x(0,0), &inc, &one, &(*rec_)(0,0), &inc); 
#elif defined PGI
  lapack::dgemv(transa, rows, cols, one, &A(0,0), lda, 
		&x(0,0), inc, one, &(*rec_)(0,0), inc); 
#elif defined  GCC
   lapack::dgemv_(&transa, &rows, &cols, &one, &A(0,0), &lda, 
 		 &x(0,0), &inc, &one, &(*rec_)(0,0), &inc); 
#endif
}


// A dgemv wrapper: Performs y = A.x + y, where x and y are vectors
void Blas_matrix::dgemm(char transa, char transb, Blas_matrix &A, Blas_matrix &B)
{
#ifdef SAFE
#endif  
  
  double one = 1.;  
  
  

  int m = this->rows();
  int n = this->cols();
  int k;
  if (transa == 'N') k = A.cols();
  else k = A.rows();
  int lda, ldb;
  if (transa == 'N') lda = m;
  else lda = k;
  if (transb == 'N') ldb = k;
  else ldb = n;
  int ldc = m;

  if(!m || !n || !k) return;
#if defined  INTEL
   lapack::dgemm(&transa, &transb, &m, &n,  &k, &one, &A(0,0), &lda, 
 		&B(0,0), &ldb, &one, &(*rec_)(0,0), &ldc); 
#elif defined PGI
  lapack::dgemm(transa, transb, m, n,  k, one, &A(0,0), lda, 
		&B(0,0), ldb, one, &(*rec_)(0,0), ldc); 
#elif defined  GCC
   lapack::dgemm_(&transa, &transb, &m, &n,  &k, &one, &A(0,0), &lda, 
 		&B(0,0), &ldb, &one, &(*rec_)(0,0), &ldc); 
#endif
}





// A dsymv wrapper: Performs y = A.x + y, A is symmetric
void Blas_matrix::dsymv(Blas_matrix &A, Blas_matrix &x)
{
#ifdef SAFE
  if((cols() != 1) || (x.cols() != 1))
    throw std::string("Blas_matrix::dsymv(): Term not a vector.");
#endif  

  double one = 1.;  
  int rows = A.rows(), inc = 1, lda = rows;
  char uplo = 'L';
#if defined  INTEL 
   lapack::dsymv(&uplo, &rows, &one, &A(0,0), &lda, 
 		 &x(0,0), &inc, &one,  &(*rec_)(0,0), &inc);
#elif defined PGI
  lapack::dsymv(uplo, rows, one, &A(0,0), lda, 
		&x(0,0), inc, one,  &(*rec_)(0,0), inc);
#elif defined  GCC 
   lapack::dsymv_(&uplo, &rows, &one, &A(0,0), &lda, 
 		 &x(0,0), &inc, &one,  &(*rec_)(0,0), &inc);
#endif  
}


// Gets the eigenvalues of the local symmetric matrix
void Blas_matrix::dsyev(char jobz, Blas_matrix& eigen_values)
{
#ifdef SAFE
  if (!size())
    throw std::string("Blas_matrix::dsyev: Cannot diagonalize 0 size matrix.");
  if (cols() != rows())
    throw std::string("Blas_matrix::dsyev: Matrix not square.");
  if ((jobz != 'N') && (jobz != 'V'))
    throw std::string("Blas_matrix::dsyev: Wrong input parameter jobz.");    
#endif
  
  eigen_values.allocate(rows());

  int n = rows(), lda = rows(), lwork = 3 * rows(), info = 0;
  double *work = new double[lwork];
  char uplo = 'L';
#if defined  INTEL  
   lapack::dsyev(&jobz, &uplo, &n, &(*rec_)(0,0), &lda, 
 		 &eigen_values(0,0), work, &lwork, &info);
#elif defined PGI
  lapack::dsyev(jobz, uplo, n, &(*rec_)(0,0), lda, 
		 &eigen_values(0,0), &info);
#elif defined  GCC  
   lapack::dsyev_(&jobz, &uplo, &n, &(*rec_)(0,0), &lda, 
 		 &eigen_values(0,0), work, &lwork, &info);
#endif  
  delete [] work;

#ifdef SAFE
  if (info)
    throw std::string("Blas_matrix::dsyev: dsyev failed.");  
#endif  

}



// Required by the Lanczos template: a linear combination of vectors
 void Blas_matrix::LinCombo(double alpha, Blas_matrix &bx,
			    double beta, Blas_matrix &by)
{
  
#ifdef SAFE
  if (bx.size() != by.size())
    throw std::string("Blas_matrix::LinCombo: Different vector sizes.");
#endif  
  
  double *first_coef, *second_coef;
  Blas_matrix *second_term;
  
  if ( (&bx == this) && (&by == this) ) {
    this->dscal(alpha + beta);
    return;
  } else if ( &bx == this ) {
    first_coef = &alpha;
    second_coef = &beta;    
    second_term = &by;
  } else if ( &by == this ) {
    first_coef = &beta;
    second_coef = &alpha;
    second_term = &bx;
  } else  {
    first_coef = &alpha;
    second_coef = &beta;    
    second_term = &by;
    this->dcopy(bx);
  }
  
  this->dscal(*first_coef);
  this->daxpy(*second_coef, *second_term);
}


void Blas_matrix::inverse()
{
#ifdef SAFE
  if (!size())
    throw std::string("Blas_matrix::inverse: Cannot invert 0 size matrix.");
  if (rows() != cols())
    throw std::string("Blas_matrix::inverse: Not a square matrix.");

#endif
  
  int m = rows();
  int n = cols();
  int lda = m;
  int info;

  int* ipiv = new int[(m > n) ? m : n];
#if defined  INTEL
  lapack::dgetrf(&m, &n, &(*rec_)(0,0), &lda, ipiv, &info);
#elif defined PGI
  lapack::dgetrf(m, n, &(*rec_)(0,0), lda, ipiv, &info);
#elif defined  GCC
  lapack::dgetrf_(&m, &n, &(*rec_)(0,0), &lda, ipiv, &info);
#endif  
#ifdef SAFE
  if (info)
    throw std::string("Blas_matrix::inverse: dgetrf (LU factorization) failed.");  
#endif  
#if defined  INTEL
  int ispec = 1;
  int n1,n2,n3,n4;

  int lwork = n * lapack::ilaenv(&ispec, "dgetri", "", &n1, &n2, &n3, &n4);

  double* work = new double[lwork];

  lapack::dgetri(&n, &(*rec_)(0,0), &lda, ipiv, work, &lwork, &info);
  delete [] work;
#elif defined PGI
  lapack::dgetri(n, &(*rec_)(0,0), lda, ipiv, &info);
#elif defined  GCC
  int ispec = 1;
  int n1,n2,n3,n4;
  char name[] = "dgetri";
  char opt[]  = "";
	
  int lwork = n * lapack::ilaenv_(&ispec, name, opt, &n1, &n2, &n3, &n4);

  double* work = new double[lwork];

  lapack::dgetri_(&n, &(*rec_)(0,0), &lda, ipiv, work, &lwork, &info);
  delete [] work;
#endif  
#ifdef SAFE
  if (info)
    throw std::string("Blas_matrix::inverse: dgetri failed.");  
#endif  

  delete [] ipiv;
}



//Submatrix methods
// Similar to above

Submatrix& Submatrix::operator=(Blas_matrix &m)
{
  
#ifdef SAFE
  if((row_slice != m.rows()) || ( col_slice != m.cols()))
    throw std::string("Submatrix::operator=(): Different matrix dimensions.");
#endif
  int inc = 1;
  int size = row_slice;
  for (unsigned col = 0; col < col_slice; col++)
#if defined  INTEL
     lapack::dcopy(&size, &m(0, col), &inc, 
 		&matrix(row_start, col + col_start), &inc); 
#elif defined PGI
    lapack::dcopy(size, &m(0, col), inc, 
		&matrix(row_start, col + col_start), inc); 
#elif defined  GCC
     lapack::dcopy_(&size, &m(0, col), &inc, 
 		&matrix(row_start, col + col_start), &inc); 
#endif  
  return *this;
}





Submatrix& Submatrix::operator=(double alpha)
{
  
  int inc = 1, no_inc = 0;
  int size = row_slice;
  for (unsigned int col = 0; col < col_slice; col++)
#if defined  INTEL
     lapack::dcopy(&size, &alpha, &no_inc, 
 		&matrix(row_start, col + col_start), &inc); 
#elif defined PGI
    lapack::dcopy(size, &alpha, no_inc, 
		&matrix(row_start, col + col_start), inc); 
#elif defined  GCC
     lapack::dcopy_(&size, &alpha, &no_inc, 
 		&matrix(row_start, col + col_start), &inc); 
#endif  
  return *this;
}


Submatrix& Submatrix::operator+=(double alpha)
{
  double one = 1.;
  int inc = 1, no_inc = 0;
  int size = row_slice;
  for (unsigned int col = 0; col < col_slice; col++)
#if defined  INTEL
     lapack::daxpy(&size, &one,&alpha, &no_inc, 
 		&matrix(row_start, col + col_start), &inc); 
#elif defined PGI
  lapack::daxpy(size, one,&alpha, no_inc, 
		&matrix(row_start, col + col_start), inc); 
#elif defined  GCC
     lapack::daxpy_(&size, &one,&alpha, &no_inc, 
 		&matrix(row_start, col + col_start), &inc); 
#endif  
  return *this;
  
}


// A dgemv wrapper: Performs y = A.x + y, where x and y are vector addresses
void Submatrix::dgemv(char transa, Blas_matrix &A, const Submatrix &x) 
{
#ifdef SAFE
  if((col_slice != 1) || (x.col_slice != 1))
    throw std::string("Submatrix::dgemv(): Matrix slice not a vector.");
  if(
     ((transa == 'N') && ((row_slice != A.rows()) || 
			  (A.cols() != x.row_slice)))
     ||
     ((transa == 'T') && ((row_slice != A.cols()) || 
			  (A.rows() != x.row_slice)))
     ) {
    throw std::string("Submatrix::dgemv(): Different matrix sizes.");
  }
#endif
  double one = 1.;
  int inc = 1, rows = A.rows(), cols = A.cols(), lda = rows;
#if defined  INTEL  
   lapack::dgemv(&transa, &rows, &cols, &one, &A(0, 0), &lda, 
 	      &x.matrix(x.row_start, x.col_start),
 		 &inc, &one, &matrix(row_start, col_start), &inc);
#elif defined PGI
  lapack::dgemv(transa, rows, cols, one, &A(0, 0), lda, 
	      &x.matrix(x.row_start, x.col_start),
		inc, one, &matrix(row_start, col_start), inc);
#elif defined  GCC  
   lapack::dgemv_(&transa, &rows, &cols, &one, &A(0, 0), &lda, 
 	      &x.matrix(x.row_start, x.col_start),
 		 &inc, &one, &matrix(row_start, col_start), &inc);
#endif  
}



void Submatrix::daxpy(double alpha, Blas_matrix& bm)
{

#ifdef SAFE  
  if( col_slice * row_slice != bm.size() )
    throw std::string("Submatrix::daxpy: Different matrix sizes.");
#endif

  int inc = 1;
  int size = row_slice;
  for (unsigned int col = 0; col < col_slice; col++)
#if defined  INTEL
     lapack::daxpy(&size, &alpha, &bm(0, col), &inc, 
 		&matrix(row_start, col + col_start), &inc);
#elif defined PGI
    lapack::daxpy(size, alpha, &bm(0, col), inc, 
		   &matrix(row_start, col + col_start), inc);
#elif defined  GCC
     lapack::daxpy_(&size, &alpha, &bm(0, col), &inc, 
 		&matrix(row_start, col + col_start), &inc);
#endif  
}




void Submatrix::daxpy(double alpha, Submatrix bm)
{
  
#ifdef SAFE  
  if((col_slice != bm.col_slice) || (row_slice != bm.row_slice))
    throw std::string("Submatrix::daxpy: Different matrix sizes.");
#endif

  int inc = 1;
  int size = row_slice;
  for (unsigned int col = 0; col < col_slice; col++)
#if defined  INTEL
     lapack::daxpy(&size, &alpha, &bm.matrix(bm.row_start, col + bm.col_start), 
 		&inc, &matrix(row_start, col + col_start), &inc);
#elif defined PGI
    lapack::daxpy(size, alpha, &bm.matrix(bm.row_start, col + bm.col_start), 
		inc, &matrix(row_start, col + col_start), inc);
#elif defined  GCC
     lapack::daxpy_(&size, &alpha, &bm.matrix(bm.row_start, col + bm.col_start), 
 		&inc, &matrix(row_start, col + col_start), &inc);
#endif  
}







// Those are almost identical to the Blas_matrix's daxpy methods
void Submatrix::daxpy(double alpha, Daxpy_argument* b)
{
  b->get_argument(alpha, *this);
}


void Submatrix::daxpy(double alpha, Rectangular_argument &r)
{
#ifdef SAFE
  if ((row_slice != r.rows()) || (col_slice != r.cols()))
    throw std::string("Submatrix::daxpy: rec Different vector sizes.");
#endif  
  
  int inc = 1, length = row_slice;
  
  for(unsigned int col = 0; col < col_slice; col++)
#if defined  INTEL
    lapack::daxpy(&length, &alpha, &r(0,col), &inc, &matrix(row_start,col+col_start), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &r(0,col), inc, &matrix(row_start,col+col_start), inc);
#elif defined  GCC
    lapack::daxpy_(&length, &alpha, &r(0,col), &inc, &matrix(row_start,col+col_start), &inc);
#endif  
}


void Submatrix::daxpy(double alpha, Transposed_argument &t)
{
#ifdef SAFE
  if ((row_slice != t.rows()) || (col_slice != t.cols()))
    throw std::string("Submatrix::daxpy: trans Different vector sizes.");
#endif  
  
  int inc = 1, inc_row = matrix.rows(), length = col_slice;
  
  for(unsigned int row = 0; row < row_slice; row++)
#if defined  INTEL
     lapack::daxpy(&length, &alpha, &t(row,0), &inc,
 		&matrix(row_start+row,col_start), &inc_row);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(row,0), inc,
		&matrix(row_start+row,col_start), inc_row);
#elif defined  GCC
     lapack::daxpy_(&length, &alpha, &t(row,0), &inc,
 		&matrix(row_start+row,col_start), &inc_row);
#endif  
}


void Submatrix::daxpy(double alpha, Upper_triangular_argument& t)
{
#ifdef SAFE
  if ((row_slice != t.rows()) || (col_slice != t.cols()))
    throw std::string("Submatrix::daxpy:  Different vector sizes.");
#endif  

  int inc = 1, inc_row = matrix.rows();

  for(unsigned int col = 0; col < col_slice; col++) {
    
    int length = col + 1;
#if defined  INTEL    
     lapack::daxpy(&length, &alpha, &t(0,col), &inc, 
 		&matrix(row_start,col_start+col), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(0,col), inc, 
		  &matrix(row_start,col_start+col), inc);
#elif defined  GCC    
     lapack::daxpy_(&length, &alpha, &t(0,col), &inc, 
 		&matrix(row_start,col_start+col), &inc);
#endif   
    length--;
#if defined  INTEL    
     lapack::daxpy(&length, &alpha, &t(0,col), &inc, 
 		&matrix(row_start+col,col_start), &inc_row);    
#elif defined PGI    
    lapack::daxpy(length, alpha, &t(0,col), inc, 
		&matrix(row_start+col,col_start), inc_row);    
#elif defined  GCC    
     lapack::daxpy_(&length, &alpha, &t(0,col), &inc, 
 		&matrix(row_start+col,col_start), &inc_row);    
#endif

  }
}



void Submatrix::daxpy(double alpha, Lower_triangular_argument& t)
{
#ifdef SAFE
  if ((row_slice != t.rows()) || (col_slice != t.cols()))
    throw std::string("Blas_matrix::daxpy:  Different vector sizes.");
#endif  

  int inc = 1, inc_row = matrix.rows();
  
  for(unsigned int row = 0; row < row_slice; row++) {
    
    int length = row_slice - row;
#if defined  INTEL    
     lapack::daxpy(&length, &alpha, &t(row,row), &inc,
		&matrix(row_start+row,col_start+row), &inc);
#elif defined PGI
    lapack::daxpy(length, alpha, &t(row,row), inc,
		&matrix(row_start+row,col_start+row), inc);
#elif defined  GCC    
     lapack::daxpy_(&length, &alpha, &t(row,row), &inc,
		&matrix(row_start+row,col_start+row), &inc);
#endif
    length--;

    if(length)
#if defined  INTEL
       lapack::daxpy(&length, &alpha, &t(row+1,row), &inc, 
 		  &matrix(row_start+row,col_start+row+1), &inc_row);    
#elif defined PGI
      lapack::daxpy(length, alpha, &t(row+1,row), inc, 
		  &matrix(row_start+row,col_start+row+1), inc_row);    
#elif defined  GCC
       lapack::daxpy_(&length, &alpha, &t(row+1,row), &inc, 
 		  &matrix(row_start+row,col_start+row+1), &inc_row);    
#endif
  }

}


void Submatrix::add_diag(Blas_matrix &x)
{
#ifdef SAFE
  if (row_slice != col_slice)
    throw std::string("Submatrix::add_diag: Not a square submatrix.");
  if (row_slice != x.size())
    throw std::string("Submatrix::add_diag: Vector length different from the diagonal length.");  
#endif

  double alpha = 1.;
  int inc = 1, diag_inc = matrix.rows() + 1;
  int size = row_slice;
#if defined  INTEL
   lapack::daxpy(&size, &alpha, &x(0, 0), &inc, 
 	      &matrix(row_start, col_start), &diag_inc);
#elif defined PGI
  lapack::daxpy(size, alpha, &x(0, 0), inc, 
	      &matrix(row_start, col_start), diag_inc);
#elif defined  GCC
   lapack::daxpy_(&size, &alpha, &x(0, 0), &inc, 
 	      &matrix(row_start, col_start), &diag_inc);
#endif  
}

void Submatrix::add_diag(double alpha)
{
#ifdef SAFE
  if (row_slice != col_slice)
    throw std::string("Submatrix::add_diag: Not a square submatrix.");
#endif

  double one = 1.;  
  int no_inc = 0, diag_inc = matrix.rows() + 1;
  int size= row_slice;
#if defined  INTEL
   lapack::daxpy(&size, &one, &alpha, &no_inc, 
 	      &matrix(row_start, col_start), &diag_inc);
#elif defined PGI
  lapack::daxpy(size, one, &alpha, no_inc, 
	      &matrix(row_start, col_start), diag_inc);
#elif defined  GCC
   lapack::daxpy_(&size, &one, &alpha, &no_inc, 
 	      &matrix(row_start, col_start), &diag_inc);
#endif  
}


// A ddot wrapper
double Submatrix::ddot(const Submatrix &x) 
{
#ifdef SAFE
  if((x.col_slice != 1) || (col_slice != 1))
    throw std::string("Submatrix::ddot(): Term not a vector.");
  if(x.row_slice != row_slice)
    throw std::string("Submatrix::ddot():  Different vector sizes.");
#endif
  
  int inc = 1;
  int size = x.row_slice;
#if defined  INTEL
   return lapack::ddot(&size, &matrix(row_start,col_start), &inc, 
 		    &x.matrix(x.row_start, x.col_start), &inc);
#elif defined PGI
  return lapack::ddot(size, &matrix(row_start,col_start), inc, 
		    &x.matrix(x.row_start, x.col_start), inc);
#elif defined  GCC
   return lapack::ddot_(&size, &matrix(row_start,col_start), &inc, 
 		    &x.matrix(x.row_start, x.col_start), &inc);
#endif  
}


//A ddot wrapper overloading the multiplication operator
double Submatrix::operator*(const Submatrix &bx) 
{
  return this->ddot(bx);  
}

#ifndef __BLAS_MATRIX_HPP__
#define __BLAS_MATRIX_HPP__

// This file contains the declaration of two classes: Submatrix and Blas_matrix.
// The Blas_matrix is intended to operate on its own private matrix employing some BLAS and LAPACK routines,
// while Submatrix was intended as its proxy which defines just a subarea of that matrix. 
// Unfortunately those classes have gone slightly out of control and their design
// may have to be rethought (when/if there's time). There are a couple of problems here:

// 1) Blas_matrix operates on a rectangular matrix argument, a contiguous chunk of memory, while Submatrix
// points to a submatrix of the Blas_matrix's matrix: the memory area is no longer
// contiguous (but is composed of contiguous parts). Thus one has to reimplement the methods that operate on it. 
// So Submatrix became in charge of operations beyond its proxy responsibilities, though this was not the original intention. 

// 2) Later on it became evident that the arguments that were to be fed to Blas_matrix's and
// Submatrix's daxpy could have even more complex memory structures (being submatrices of the
// integral tables). So there had to be included a set of classes: class Daxpy_argument,
// class Rectangular_argument, class Transposed_argument, class Upper_triangular_argument,
// class Lower_triangular_argument (see integral_blocks.hpp). Those required a specific way in which the
// BLAS's daxpy could be applied to them, so everything got messier.

// The current code should be bug free, however it's far from general (and elegant).
// There should be a nicer way to handle that dependence on the memory structure
// of the arguments.


#include "matrices.hpp"

class Daxpy_argument;
class Rectangular_argument;
class Transposed_argument;
class Upper_triangular_argument;
class Lower_triangular_argument;

class Blas_matrix;





class Submatrix {

  Blas_matrix &matrix;
  unsigned int row_start, row_slice;
  unsigned int col_start, col_slice;
  
  Submatrix(Blas_matrix &m, unsigned a, unsigned b, unsigned c, unsigned d)  
    : matrix (m), row_start(a), row_slice(b), col_start(c), col_slice(d){}

public:


  double operator*(const Submatrix &bx);
  double ddot(const Submatrix &x);
  
  Submatrix& operator=(Blas_matrix &m);
  Submatrix& operator=(double alpha);
  Submatrix& operator+=(double alpha);

  void daxpy(double alpha, Blas_matrix& x);
  void daxpy(double alpha, Submatrix bm);
  void daxpy(double alpha, Rectangular_argument &r);
  void daxpy(double alpha, Transposed_argument &r);
  void daxpy(double alpha, Upper_triangular_argument &r);
  void daxpy(double alpha, Lower_triangular_argument &r);
  void daxpy(double alpha, Daxpy_argument* r);



  void add_diag(Blas_matrix &x);
  void add_diag(double alpha);

  void dgemv(char transa, Blas_matrix &x, const Submatrix &y);

  friend class Blas_matrix;  
};










class Blas_matrix {
  
  Rectangular_matrix<double>* rec_;

public:
  
  Blas_matrix(unsigned int rows = 0, unsigned int cols = 1, double *mat = 0);

  ~Blas_matrix() { delete rec_;}
  inline unsigned int rows() const {return rec_->rows();}
  inline unsigned int cols() const {return rec_->cols();}
  inline unsigned int size() const {return rec_->size();}

  inline void allocate(unsigned int rows, unsigned int cols = 1)
  { delete rec_;  
    rec_ = new Rectangular_matrix<double>(rows, cols); }
  
  // Required by the Lanzcos template:
  Blas_matrix(const Blas_matrix &m);
  double operator||(unsigned int l);
  double operator*(Blas_matrix &m);
  Blas_matrix& operator/=(double alpha);
  void LinCombo(double alpha, Blas_matrix &x, double beta, Blas_matrix  &y);

  
  //indexing
  inline double& operator()(unsigned int row, unsigned int col = 0)
  { return (*rec_)(row, col); }
  
  inline double& operator()(unsigned int row, unsigned int col = 0) const
  { return (*rec_)(row, col); } // Added just because of the copy constuctor 

  inline Submatrix operator()(unsigned row_start, unsigned row_slice,
			      unsigned col_start, unsigned col_slice)
  { return Submatrix(*this, row_start, row_slice, col_start, col_slice); }

  //dcopy wrappers
  void dcopy(Blas_matrix &m);
  void dcopy(Submatrix &m);
  void dcopy(double alpha);
  Blas_matrix& operator=(Blas_matrix &m);
  Blas_matrix& operator=(double alpha);

  //ddot wrappers
  double ddot(Blas_matrix &x);
  double ddot(const Submatrix &x);
  double operator*(const Submatrix &x);
  
  //dscal wrappers
  void dscal(double alpha);
  Blas_matrix& operator*=(double alpha);
  
  //a daxpy wrapper
  void daxpy(double alpha, Blas_matrix &x);
  
  //with submatrix arguments
  void daxpy(double alpha, Daxpy_argument* b);
  void daxpy(double alpha, Rectangular_argument& r);
  void daxpy(double alpha, Transposed_argument& t);
  void daxpy(double alpha, Upper_triangular_argument& t);
  void daxpy(double alpha, Lower_triangular_argument& t);

  void add_diag(Blas_matrix &x);
  void add_diag(double alpha);


  //dgemv wrappers
  void dgemv(char transa, Blas_matrix &A, Blas_matrix &x);
  void dgemm(char transa, char transx, Blas_matrix &A, Blas_matrix &x);
  void dsymv(Blas_matrix &A, Blas_matrix &x);
  
  // lapack diagonalizer
  void dsyev(char jobz, Blas_matrix& eigen_values);

  

  //dgetri
  void inverse();

  // misc
  void print();
};




#endif //#ifndef __BLAS_MATRIX_HPP__

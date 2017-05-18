#ifndef __INTEGRAL_BLOCKS_HPP__
#define __INTEGRAL_BLOCKS_HPP__

// This file declares a set of proxy classes that reference areas of the two-electron integral
// tables (see integral_table.hpp) corresponding to the integral blocks Vij,k
// Aij and Bij as defined in Chem.Phys. 329, p.13. These blocks are to be fed to
// the BLAS daxpy routine, although they are not generally monolith contiguous chunks of memory. 
// The purpose of the classes here is to define the memory "structure" of the blocks
// in order for Blas_matrix (see blas_matrix.hpp) to know how to apply the BLAS daxpy routine.


#include "matrices.hpp"
#include "blas_matrix.hpp"



// The abstract base class of all integral block classes.
// It declares the virtual get_argument methods which are used
// by Blas_matrix to determine which type of integral block
// is the current argument of daxpy (see blas_matrix.hpp)
class Daxpy_argument: public Matrix<double> {
public:
  virtual void get_argument(double alpha, Blas_matrix& b) = 0;
  virtual void get_argument(double alpha, Submatrix& b) = 0;
};


// An abstract class that defines rectangular integral blocks whose
// columns are contiguous chunks 
class Rectangular_argument: public Daxpy_argument {
protected:

  unsigned int row_start_, row_slice_;
  unsigned int col_start_, col_slice_;
  
public:
  
  virtual void get_argument(double alpha, Blas_matrix& b) { b.daxpy(alpha, *this); }
  virtual void get_argument(double alpha, Submatrix& b)   { b.daxpy(alpha, *this); }

  Rectangular_argument(unsigned int row_start, unsigned int row_slice,
		       unsigned int col_start, unsigned int col_slice)
    : row_start_(row_start), row_slice_(row_slice),
      col_start_(col_start), col_slice_(col_slice) {}
  
  virtual unsigned int rows() {return row_slice_;}
  virtual unsigned int cols() {return col_slice_;}
  virtual unsigned int size() {return row_slice_ * col_slice_;}
};


// A rectangular integral block that is a submatrix of a larger matrix.
// Its columns are parts of the columns of the larger matrix.
// Such types of blocks are the Vij,k blocks (a vector), the Aij blocks for i<j.
class Rectangular_submatrix: public Rectangular_argument {
  
  Matrix<double>& mat_;
  
public:
  
  Rectangular_submatrix(Matrix<double>& mat, 
			unsigned int row_start, unsigned int row_slice,
			unsigned int col_start, unsigned int col_slice)
    : mat_(mat), Rectangular_argument(row_start, row_slice,
				      col_start, col_slice)
  {
    
//     if ((row_start + row_slice > mat.rows())
// 	|| (col_start + col_slice > mat.cols())) 
//     	throw std::string("Rec Submatrix out of boundaries.\n");
      
  }
  
  virtual double& operator()(unsigned int row, unsigned int col)
  {return mat_(row_start_ + row, col_start_ + col);}
  
};


// A rectangular integral block which is in fact a contiguous chunk of memory. It
// is a single column (or its part) of the integral table that must be
// treated as a rectangular matrix.
// Such are the Bij blocks for non-totally symmetric virtual orbitals where the symmetry
// of the row index is lower than the symmetry of the column index (see integral_table.cpp)
class Rectangular_proxymatrix: public Rectangular_argument {
  
  double* start_;
  
public:
  
  Rectangular_proxymatrix(Matrix<double>& mat, 
			   unsigned int row_start, unsigned int row_slice,
			   unsigned int col_start, unsigned int col_slice)
    : Rectangular_argument(row_start, row_slice,
			   col_start, col_slice)
  {
//     if (row_start + row_slice*col_slice > mat.rows()) 
      
//       throw std::string("Rec Submatrix out of boundaries.\n");
    
    start_ = &mat(row_start,col_start);
  }

  Rectangular_proxymatrix(double* mat, unsigned int row_slice,
			  unsigned int col_slice)
    : Rectangular_argument(0, row_slice,
			   0, col_slice)
  { start_ = mat; }

  
  virtual double& operator()(unsigned int row, unsigned int col)
  {return start_[row + col*row_slice_];}
  
};

// An abstract class that defines rectangular integral blocks whose
// rows are contiguous chunks of memory
class Transposed_argument: public Daxpy_argument {
  
protected:

  unsigned int row_start_, row_slice_;
  unsigned int col_start_, col_slice_;

public:
  
  virtual void get_argument(double alpha, Blas_matrix& b) {

    b.daxpy(alpha, *this);}
  virtual void get_argument(double alpha, Submatrix& b) { b.daxpy(alpha, *this);}


  Transposed_argument(unsigned int row_start, unsigned int row_slice,
		      unsigned int col_start, unsigned int col_slice)
    : row_start_(row_start), row_slice_(col_slice),
      col_start_(col_start), col_slice_(row_slice) 
  {  }
  
  virtual unsigned int rows() {return row_slice_;}
  virtual unsigned int cols() {return col_slice_;}
  virtual unsigned int size() {return row_slice_ * col_slice_;}
};


// A rectangular integral block that is a transposed submatrix of a larger matrix.
// Its rows are parts of the columns of the larger matrix.
// Such types of blocks are the Aij blocks for j>i.
class Transposed_submatrix: public Transposed_argument {
  
  Matrix<double>& mat_;

public:
  
  Transposed_submatrix(Matrix<double>& mat, 
		       unsigned int row_start, unsigned int row_slice,
		       unsigned int col_start, unsigned int col_slice)
    : mat_(mat), Transposed_argument(row_start, row_slice,
				    col_start, col_slice) 
  {
//     if ((row_start + row_slice > mat.rows())
// 	|| (col_start + col_slice > mat.cols()))

// 	throw std::string("Trans Submatrix out of boundaries.\n");
  }
    
  virtual double& operator()(unsigned int row, unsigned int col)
  {return mat_(row_start_ + col, col_start_ + row);}
  
};




// A rectangular integral block that
// references a part of a column of the integral table and must be
// treated as a transposed matrix.
// Such are the Bij blocks for non-totally symmetric virtual orbitals where the symmetry
// of the row index is bigger than the symmetry of the column index (see integral_table.cpp)
class Transposed_proxymatrix: public Transposed_argument {
  
  double* start_;

public:
  
  Transposed_proxymatrix(Matrix<double>& mat, 
			  unsigned int row_start, unsigned int row_slice,
			  unsigned int col_start,unsigned int col_slice)
    : Transposed_argument(row_start, row_slice,
			  col_start, col_slice) 
  {
//     if (row_start + row_slice*col_slice > mat.rows())
      
//       throw std::string("Trans Submatrix out of boundaries.\n");
    
    start_ = &mat(row_start,col_start);
    
  }
  
  Transposed_proxymatrix(double* mat, 
			 unsigned int row_slice,
			 unsigned int col_slice)
    : Transposed_argument(0, row_slice,
			  0, col_slice) 
  {
    start_ = mat;
  }


  virtual double& operator()(unsigned int row, unsigned int col)
  {return start_[col + row * col_slice_];}
  
};




// A class that defines an upper triangular area (lying on the diagonal)
// of a larger triangular matrix. Its columns are parts of the
// columns of the larger matrix and are contiguous chunks.
// Such are the Aii blocks.
class Upper_triangular_argument: public Daxpy_argument {

  Triangular_matrix<double>& mat_;  
  unsigned int row_start_, col_start_, row_slice_;

  static unsigned int sum(unsigned int a) {return a*(a+1)/2;}
  
public:
  
  virtual void get_argument(double alpha, Blas_matrix& b) {
    b.daxpy(alpha, *this);}
  virtual void get_argument(double alpha, Submatrix& b) { b.daxpy(alpha, *this);}


  Upper_triangular_argument(Triangular_matrix<double>& mat,
				 unsigned int row_start, unsigned int col_start, 
				 unsigned int row_slice)
    : mat_(mat), row_start_(row_start), row_slice_(row_slice),
      col_start_(col_start) 
  {
//     if (row_start != col_start)
//       throw std::string("The submatrix is not triangular.\n");
    
//     if (row_start + row_slice > mat_.rows())
//       throw std::string("The submatrix is not triangular.\n");
  }
  
  virtual unsigned int rows() {return row_slice_;}
  virtual unsigned int cols() {return row_slice_;}
  virtual unsigned int size() {return sum(row_slice_);}
 
  virtual double& operator()(unsigned int row, unsigned int col)
  {return mat_(row_start_ + row, col_start_ + col);}
};


// A triangular integral block that
// references a part of a single column of the integral table and must be
// treated as a lower triangular matrix. Its columns are contiguous chunks.
// Such are the Bij blocks for totally symmetric virtual orbitals
class Lower_triangular_argument: public Daxpy_argument {
  
  double* start_;

  unsigned int row_start_, col_start_, row_slice_;
  
  static unsigned int sum(unsigned int a) {return a*(a+1)/2;}
  static void swap(unsigned int &a, unsigned int& b) {unsigned int t=a;a=b;b=t;}
  
public:
  
  virtual void get_argument(double alpha, Blas_matrix& b) {
    b.daxpy(alpha, *this);}
  virtual void get_argument(double alpha, Submatrix& b) { b.daxpy(alpha, *this);}

  
  
  Lower_triangular_argument(Rectangular_matrix<double>& mat, 
				 unsigned int row_start, unsigned int col_start, 
				 unsigned int row_slice) :
    row_start_(row_start), row_slice_(row_slice),
    col_start_(col_start) {
//     if (row_start + sum(row_slice) > mat.rows())
//       throw std::string("The submatrix is not triangular.\n");

    start_ = &mat(row_start, col_start);
  }


  virtual unsigned int rows() {return row_slice_;}
  virtual unsigned int cols() {return row_slice_;}
  virtual unsigned int size() {return sum(row_slice_);}
  
  virtual double& operator()(unsigned int row, unsigned int col)
  {
    if (col > row) swap(row,col);
    return start_[row + size() - sum(cols() - col) - col] ;
  }
  
};



#endif //#ifndef __INTEGRAL_BLOCKS_HPP__

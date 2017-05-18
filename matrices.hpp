#ifndef __MATRICES_HPP__
#define __MATRICES_HPP__



template <class T>
class Matrix {
public:
  virtual ~Matrix() {}
  virtual unsigned int rows() = 0;
  virtual unsigned int cols() = 0;
  virtual unsigned int size() = 0;
  virtual T& operator()(unsigned int i, unsigned int j) = 0;
};



template <class T> 
class Rectangular_matrix: public Matrix<T> {
  
  unsigned int rows_, cols_;
  T* data_;
  
public:
  
  Rectangular_matrix(unsigned int rows, unsigned int cols) : rows_(rows), cols_(cols) 
  { 
    if (!rows_ || !cols_) rows_ = cols_ = 0; 
    data_ = new T[rows_ * cols_];
  }
  
  virtual ~Rectangular_matrix() {delete [] data_;}
  virtual unsigned int rows() {return rows_;}
  virtual unsigned int cols() {return cols_;}
  virtual unsigned int size() {return rows_ * cols_;}
  inline virtual T& operator()(unsigned int row, unsigned int col)
  {/*add boundery checks*/return data_[col * rows_ + row];}

};


template <class T>
class Triangular_matrix: public Matrix<T> {
  
  unsigned int rows_;
  T* data_;
  
  inline unsigned int sum(unsigned int a) {return (a*(a+1))>>1;}
  inline void swap(unsigned int& a, unsigned int& b) {unsigned int t=a;a=b;b=t;}
  
  inline int abs(int ecx)
  {
    int ebx = ecx;
    ecx = ecx >> 31;
    ebx = ebx ^ ecx;
    ebx -= ecx;
    return ebx;
  }
public:

  explicit Triangular_matrix(unsigned int rows) : rows_(rows)
  { data_ = new T[sum(rows)];}

  virtual ~Triangular_matrix() {delete [] data_;}
  
  virtual unsigned int rows() {return rows_;}
  virtual unsigned int cols() {return rows_;}
  virtual unsigned int size() {return sum(rows_);}
  inline virtual T& operator()(unsigned int row, unsigned int col)
  //{if (row > col) swap(row, col); return data_[sum(col) + row];}
  { unsigned s = row+col;
    col = (s + abs(row-col))>>1;
    row = s - col;
    return data_[sum(col) + row];}
};


#endif //#ifndef __MATRICES_HPP__

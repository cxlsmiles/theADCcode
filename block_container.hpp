#ifndef  __BLOCK_CONTAINER_HPP__
#define  __BLOCK_CONTAINER_HPP__

// This file declares two templates whose purpose is to store pointers
// to the blocks of integrals Vij,k: (rk|ij); Aij: (ri|sj); Bij: (rs|ij) (see Chem.Phys. 329, 13)
// For Vij,k and Bij,  i >= j, so the matrix that stores the pointers must be triangular.
// Aij and Bij, as defined in the article, can be further separated by symmetry.
// So a third parameter is needed (the symmetry of one of the virtual orbitals) to 
// specify the sub-block.

#include <vector>
#include "matrices.hpp"

#define CHECK
#undef CHECK

// Declaration:
template <class T> 
class Triangular_container {
  
  std::vector< Triangular_matrix<T*>* >  block_; 
    
public:
  
  Triangular_container(unsigned num_blocks, unsigned row_size);
  ~Triangular_container();
  T*& operator()(unsigned a, unsigned b, unsigned c);
};

template<class T> 
class Rectangular_container {
  
  std::vector< Rectangular_matrix<T*>* >  block_;

public:
  
  Rectangular_container(unsigned num_blocks, unsigned row_size);
  ~Rectangular_container();
  T*& operator()(unsigned a, unsigned b, unsigned c);
};


// Definition:
template <class T> 
Triangular_container<T>::Triangular_container(unsigned num_blocks, unsigned row_size) 
{
  
  for(unsigned i = 0; i < num_blocks; i++)
    block_.push_back(new Triangular_matrix<T*>(row_size));
  
  for(unsigned i = 0; i < block_.size(); i++)
    for(unsigned col = 0; col < block_[i]->cols(); col++)
      for(unsigned row = 0; row <= col; row++)
	(*block_[i])(row,col) = 0;
}

template<class T> 
Triangular_container<T>::~Triangular_container() 
{
  for(unsigned i = 0; i < block_.size(); i++) {
    for(unsigned col = 0; col < block_[i]->cols(); col++)
      for(unsigned row = 0; row <= col; row++)
	delete (*block_[i])(row,col);
    
    delete block_[i];
  }
}

template<class T> 
T*& Triangular_container<T>::operator()(unsigned a, unsigned b, unsigned c)
{
#ifdef CHECK  
  if ((block_.size() <= c) || (block_[c]->rows() <= a) 
      || (block_[c]->cols() <= b)) 
    throw std::string("Triangular_container::operator(): index out of bounds.\n");
#endif

  return block_[c]->operator()(a,b);
}

template<class T> 
Rectangular_container<T>::Rectangular_container(unsigned num_blocks, unsigned row_size)
{
  
  for(unsigned i = 0; i < num_blocks; i++)
    block_.push_back(new Rectangular_matrix<T*>(row_size, row_size));

  
  for(unsigned i = 0; i < block_.size(); i++)
    for(unsigned col = 0; col < block_[i]->cols(); col++)
      for(unsigned row = 0; row < block_[i]->rows(); row++)
	(*block_[i])(row,col) = 0;
}

template<class T> 
Rectangular_container<T>::~Rectangular_container()
{
  
  for(unsigned i = 0; i < block_.size(); i++) {
    for(unsigned col = 0; col < block_[i]->cols(); col++)
      for(unsigned row = 0; row < block_[i]->rows(); row++)
	delete (*block_[i])(row,col);
    delete block_[i];
  }
}

template<class T> 
T*& Rectangular_container<T>::operator()(unsigned a, unsigned b, unsigned c)
{
#ifdef CHECK  
  if ((block_.size() <= c) || (block_[c]->rows() <= a) 
      || (block_[c]->cols() <= b)) 
    throw std::string("Rectangular_container::operator(): index out of bounds.\n");
#endif

  return block_[c]->operator()(a,b);  
}

#endif //#ifndef  __BLOCK_CONTAINER_HPP__


#include "adc2_matrix.hpp"
#include "triplet.hpp"
#include "singlet.hpp"
#include "analysis/adc_analyzer.hpp"
#include "../blas_matrix.hpp"
#include "config.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;


// The file contains the implementations of the member functions of the Adc2_matrix class.

// There are two important functionalities here: building the whole matrix (mostly for testing purposes)
// and performing matrix vector multiplication (with a set of vectors) without building the whole matrix.
// The latter is used by the Lanczos diagonalizer.

// In both cases the matrix relies on the member blocks_ which supplies the building blocks. Those can either be
// computed for singlet or triplet (see singlet.hpp and triplet.hpp). Thus the implementation is defined for spin 1 or 3 only.
// 

// The matrix-vector multiply step (performed here) is the highest level at which parallelization 
// can be applied since the Lanczos procedure is iterative and each new step depends
// on the results of the preceding.  For shared memory parallelization, 
// the iterations over the blocks of the matrix had to be changed to avoid
// competition between the threads. Since the matrix is symmetric, most blocks are used twice 
// in the matrix vector multiplication, see Fig.2 in Chem.Phys. 329, 11. Therefore, 
// different threads may compete for accessing overlapping areas in the resulting vector.

// The loop order was changed to ensure that won't happen. This was done by iterating
// over the anti-diagonals of the matrix. The cells of each anti-diagonal can be processed parallelly without competition
// since those reside always in different columns/rows even when transposed.
// For a different solution see Tarantelli's article, page 16,
// however, the suggested parallelization there is in a lower implementation level:  in the computation
// of the singlet/triplet blocks.


Adc2_matrix::Adc2_matrix(SCF_data_reader& phis, Integral_table& tab, unsigned int sym, unsigned int spin) :
  phis_(phis), symmetry_(sym), begin_ij_(0), begin_jii_(0), begin_ijk_(0), count_mult(0)
{

  switch (spin) {
    //Singlet
  case 0:
    multiplet_functions_ = 2;
    blocks_ = new Singlet(phis_, tab);
    break;
    //Triplet
  case 2:
    multiplet_functions_ = 3;
    blocks_ = new Triplet(phis_, tab);
    break;
  default:
    throw string("Only singlets and triplets for ADC2 double ionization\n");
  }

  // Determine the configurations for the given symmetry and spin.
  add_ii_configs();
  add_ij_configs();
  add_jiir_configs();
  add_ijkr_configs();
  
  //print_configs();
 
}

int Adc2_matrix::accept_analyzer(ADC_analyzer& analyst)
{
  return analyst.analyze(*this);
}

Adc2_matrix::~Adc2_matrix()
{
  delete blocks_;
}


// Adds closed shell 2h configurations, |ii>
void Adc2_matrix::add_ii_configs()
{

  //no closed shell configurations for triplet
  if(multiplet_functions_ == 3) return;
  //|ii> can be only totally symmetric
  if(symmetry_) return;
  
  for(unsigned int i = 0; i < phis_.number_occupied(); i++) {
    
    Config c(i, i);
    configs_.push_back(c);
  }
  
  begin_ij_ =  configs_.size();
}

// Adds open shell 2h configurations, |ij>
void Adc2_matrix::add_ij_configs()
{
  for(unsigned int i = 0; i < phis_.number_occupied(); i++)
    for(unsigned int j = 0; j < i; j++) {
      
      if (symmetry_ != phis_.irrep_product(phis_.irrep(i), phis_.irrep(j)))
	continue;
      
      Config c(i, j);
      configs_.push_back(c);
    }
  
  begin_jii_ =  configs_.size();
  
}

//Adds 2h 1p configurations, type I, |jiir>
void Adc2_matrix::add_jiir_configs()
{
  
  jii_.push_back(begin_jii_);
  for(unsigned int j = 0; j < phis_.number_occupied(); j++) 
    for(unsigned int i = 0; i < phis_.number_occupied(); i++) {
      
      if (i == j)
	continue;
      
      bool jii_exists = false;
      for(unsigned int r = phis_.number_occupied(); r < phis_.number_orbitals(); r++) {
	
	if(symmetry_ != phis_.irrep_product(phis_.irrep(j),phis_.irrep(r)))
	  continue;
	
	Config c(j, i, i, r);
	configs_.push_back(c);
	
	jii_exists = true;
      }

      if (jii_exists) jii_.push_back(configs_.size());
    }

  jii_.pop_back();
  begin_ijk_ = configs_.size();
}

//Adds 2h 1p configurations, type II, |ijkr>     
void Adc2_matrix::add_ijkr_configs()
{
  ijk_.push_back(begin_ijk_);
  for(unsigned int i = 0; i < phis_.number_occupied(); i++) 
    
    for(unsigned int j = 0; j < i; j++) 
            
      for(unsigned int k = 0; k < j; k++) {

	bool ijk_exists = false;
	for (unsigned int type = 0;
	     type < multiplet_functions_; type++) 
	  
	  for(unsigned int r = phis_.number_occupied(); r < phis_.number_orbitals(); r++) {
	    
	    if (symmetry_ != 
		phis_.irrep_product(phis_.irrep_product(phis_.irrep_product(phis_.irrep(i),phis_.irrep(j)),phis_.irrep(k)),phis_.irrep(r)))
	      continue;
	    
	    ijk_exists = true;
	    Config c(i, j, k, r, type);
	    configs_.push_back(c);
	  }
	if (ijk_exists) {
	  
	  ijk_.push_back(configs_.size());
	  
	}
      }
  
  ijk_.pop_back();
  
}


// Builds the whole matrix
void Adc2_matrix::build_matrix(Blas_matrix &adc2_matrix)
{

  adc2_matrix.allocate(size(), size());
  adc2_matrix = 0.;

  build_submatrix_ii_jj(adc2_matrix);
  build_submatrix_ij_kk(adc2_matrix);
  build_submatrix_ij_kl(adc2_matrix);
  build_submatrix_lkk_ii(adc2_matrix);
  build_submatrix_lkk_ij(adc2_matrix);
  build_submatrix_klm_ii(adc2_matrix);
  build_submatrix_klm_ij(adc2_matrix);
  build_submatrix_jii_lkk(adc2_matrix);
  build_submatrix_ijk_mll(adc2_matrix);
  build_submatrix_ijk_lmn(adc2_matrix);

}


// Performs the matrix vector multiply step. This member function is defined
// according to the requirements of the Lanczos diagonalizer.
// The parameter count_mult counts how many times the routine has been called.
// It multiplies the main/main and the sat/main parts of the matrix with the input vectors
// only during the first two calls/Lanczos iterations. This ensures that a subspace iteration Lanczos is performed.
// For the benefits of that see Tarantelli's article, section 3.3.
// Note that the initial guess of the Lanczos procedure *must* be a set of Cartesian vectors
// with size equal to the main part of the matrix.
int Adc2_matrix::operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count)
{  


  if (count_mult < 2*begin_jii_) {

    multiply_submatrix_ii_jj(block_out, block_in, count);
    multiply_submatrix_ij_kk(block_out, block_in, count);
    multiply_submatrix_ij_kl(block_out, block_in, count);

    multiply_submatrix_lkk_ii(block_out, block_in, count);
    multiply_submatrix_lkk_ij(block_out, block_in, count);
    multiply_submatrix_klm_ii(block_out, block_in, count);
    multiply_submatrix_klm_ij(block_out, block_in, count);
  }
  
  multiply_submatrix_jii_lkk(block_out, block_in, count);
  multiply_submatrix_ijk_mll(block_out, block_in, count);
  multiply_submatrix_ijk_lmn(block_out, block_in, count);

  
  count_mult += count;

  return count;

}


// Builds the <2h, closed shell|H|2h, closed shell> block
void Adc2_matrix::build_submatrix_ii_jj(Blas_matrix &adc2_matrix)
{
  
  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;
  
  double element;
  
  for(unsigned int row = 0; row < begin_ij_; row++) 
    for(unsigned int col = 0; col < begin_ij_; col++) {
      
      blocks_->block_ii_jj(configs_[row], configs_[col], element);
      adc2_matrix(row, col) = element;
      
    }
}

// Multiplies the <2h, closed shell|H|2h, closed shell> block with the relevant part of the input vectors.
inline void 
Adc2_matrix::multiply_submatrix_ii_jj(Blas_matrix** block_out, 
				      Blas_matrix** block_in, int count)
{
  
  if(multiplet_functions_ == 3) return;
  if (symmetry_) return;
  

  double element;

  for(unsigned int row = 0; row < begin_ij_; row++) 
    for(unsigned int col = 0; col <= row ; col++) {
      
      blocks_->block_ii_jj(configs_[row], configs_[col], element);
      
      Blas_matrix **vec_in = block_in, **vec_out = block_out;
      
      for(int vec = 0; vec < count; vec++) {
	(**vec_out)(row) += element * (**vec_in)(col);

	if (row != col)
	  (**vec_out)(col) += element * (**vec_in)(row);

	vec_in++; vec_out++; 

      }

    }
}


//Builds the <2h, open shell|H|2h, closed shell> block
void Adc2_matrix::build_submatrix_ij_kk(Blas_matrix &adc2_matrix)
{
  
  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;
  
  double element;

  for(unsigned int row = begin_ij_; row < begin_jii_; row++) 
    for(unsigned int col = 0; col < begin_ij_; col++) {
      blocks_->block_ij_kk(configs_[row], configs_[col], element);
      adc2_matrix(row, col) = adc2_matrix(col, row) = element;
      
    }
}

// Multiplies the <2h, open shell|H|2h, closed shell> block  with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_ij_kk(Blas_matrix** block_out, 
				       Blas_matrix** block_in, int count)
{
  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;

  double element;

  for(unsigned row = begin_ij_; row < begin_jii_; row++) 
    for(unsigned col = 0; col < begin_ij_ ; col++) {
      
      blocks_->block_ij_kk(configs_[row], configs_[col], element);
      Blas_matrix **vec_in = block_in, **vec_out = block_out;
      
      for(int vec = 0; vec < count; vec++) {
	(**vec_out)(row) += element * (**vec_in)(col);
	if (row != col)
	  (**vec_out)(col) += element * (**vec_in)(row);
	vec_out++, vec_in++;
      }
    }
}

// Builds the <2h, open shell|H|2h, open shell> block
void Adc2_matrix::build_submatrix_ij_kl(Blas_matrix &adc2_matrix)
{
  double element;
  
  for(unsigned int row = begin_ij_; row < begin_jii_; row++) 
    for(unsigned int col = begin_ij_; col < begin_jii_; col++) {
      
      blocks_->block_ij_kl(configs_[row], configs_[col], element);
      adc2_matrix(row, col) = element;
      
    }
}


// Multiplies the <2h, open shell|H|2h, open shell> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_ij_kl(Blas_matrix** block_out, 
					 Blas_matrix** block_in, int count)
{
  double element;
  
  for(unsigned row = begin_ij_; row < begin_jii_; row++) {
    for(unsigned col = begin_ij_; col <= row ; col++) {
      
      blocks_->block_ij_kl(configs_[row], configs_[col], element);
      
      Blas_matrix **vec_in = block_in, **vec_out = block_out;
      for(int vec = 0; vec < count; vec++) {
	(**vec_out)(row) += element * (**vec_in)(col);
	if (row != col)
	  (**vec_out)(col) += element * (**vec_in)(row);
	vec_out++; vec_in++;
      }
    }
  }
}

// Builds the <3h1p,I|H|2h, closed shell> block
void Adc2_matrix::build_submatrix_lkk_ii(Blas_matrix &adc2_matrix)
{

  if(multiplet_functions_ == 3) return;
  if (symmetry_) return;

  Blas_matrix block;
  
  for(unsigned int row_inc, row = begin_jii_; 
      row < begin_ijk_; row += row_inc) {
    
    for(unsigned int col = 0; col < begin_ij_; col++) 
      
      if (blocks_->block_lkk_ii(configs_[row], configs_[col], block))
	
	adc2_matrix(row, block.rows(), col, 1) =  block;
    
    row_inc =  blocks_->size_vir_group(configs_[row].vir);   
  }  
}


// Multiplies the <3h1p,I|H|2h, closed shell> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_lkk_ii(Blas_matrix** block_out, 
					  Blas_matrix** block_in, int count)
{

  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;

  unsigned int i, j;

  unsigned int row_dim  =  jii_.size();
  unsigned int col_dim =  begin_ij_;
  
  // Return if one of the dimensions doesn't exist
  if (!row_dim || !col_dim) return;
  row_dim--; col_dim--;

#pragma omp parallel private(i,j)  
  {
    
    Blas_matrix block;
    // Loop over the anti-diagonals
    for (i= 0; i <= row_dim + col_dim; i++ ) {
      
      // Find the row and column index of the starting element (lying on the diagonal)
      // and the size of the anti-diagonal
      unsigned int row_ind; 
      unsigned int adiag_size;
      unsigned int lower_dim = row_dim < col_dim ? row_dim : col_dim;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = row_dim + col_dim - i < lower_dim ?   
	  row_dim + col_dim - i + 1  : lower_dim + 1; 
      }else {
	row_ind = i;
	adiag_size = i < lower_dim ? i + 1 : lower_dim + 1; 
      }
      unsigned int col_ind  = i - row_ind;
      
      // Loop from the middle of the anti-diagonal downwards
#pragma omp  for  schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     
	
	unsigned int row = jii_[row_ind - j];
	unsigned int col = col_ind + j;
	
	if (blocks_->block_lkk_ii(configs_[row], configs_[col], block)) { 
	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;

	  unsigned int block_rows =  block.rows();
	  
	  for(int vec = 0; vec < count; vec++) {
	    
	    (**vec_out)(row, block_rows, 0, 1).daxpy((**vec_in)(col),block);
	    
	    (**vec_out)(col) +=  
	      block
	      * (**vec_in)(row, block_rows, 0, 1);
	    
	    vec_out++; vec_in++;
	  }
	}
      }
    }    
  }
}


// Builds the <3h1p,I|H|2h, open shell> block
void Adc2_matrix::build_submatrix_lkk_ij(Blas_matrix &adc2_matrix)
{

  Blas_matrix block;
  
  for(unsigned int row_inc, row = begin_jii_; 
      row < begin_ijk_; row += row_inc) {
    
    for(unsigned int col = begin_ij_; col < begin_jii_; col++) 
      
      if (blocks_->block_lkk_ij(configs_[row], configs_[col], block)) 
	
	adc2_matrix(row, block.rows(), col, 1) =  block;
    
    
    row_inc =  blocks_->size_vir_group(configs_[row].vir);
  }  
}


// Multiplies the <3h1p,I|H|2h, open shell> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_lkk_ij(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{
  unsigned int i,j;

  unsigned int row_dim  =  jii_.size();
  unsigned int col_dim = begin_jii_ - begin_ij_;

  if (!row_dim || !col_dim) return;
  row_dim--; col_dim--;

#pragma omp parallel private(i,j)  
  {
    Blas_matrix block;
    
    for (i = 0; i <= row_dim + col_dim; i++ ) {
   
      unsigned int row_ind; 
      unsigned int adiag_size;
      unsigned int lower_dim = row_dim < col_dim ? row_dim : col_dim;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = row_dim + col_dim - i < lower_dim ?   row_dim + col_dim - i + 1  : lower_dim + 1; 
      }else {
	row_ind = i;
	adiag_size = i < lower_dim ? i + 1 : lower_dim + 1; 
      }
      unsigned int col_ind  = i - row_ind;

#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     
	
	unsigned int row = jii_[row_ind - j];
	unsigned int col =  begin_ij_ + col_ind + j;
	
	
	if (blocks_->block_lkk_ij(configs_[row], configs_[col], block)) { 
	  
	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  
	  unsigned int block_rows =  block.rows();
	  
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0, 1).daxpy((**vec_in)(col), block);
	    (**vec_out)(col) +=  
	      block
	      * (**vec_in)(row, block_rows, 0, 1);
	    vec_out++; vec_in++;
	  }
	}
      }
      
    }
  }
}

// Builds the <3h1p,II|H|2h, closed shell>
void Adc2_matrix::build_submatrix_klm_ii(Blas_matrix &adc2_matrix)
{

  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;

  Blas_matrix block;
  
  for(unsigned int row_inc, row = begin_ijk_; 
      row < configs_.size(); row += row_inc) {
    
    for(unsigned int col = 0; col < begin_ij_; col++) 
      
      if (blocks_->block_klm_ii(configs_[row], configs_[col], block))
	
	adc2_matrix(row, block.rows(), col , 1) = block;
    
    row_inc = multiplet_functions_
      *  blocks_->size_vir_group(configs_[row].vir);
    
  }
}

// Multiplies the <3h1p,II|H|2h, closed shell> with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_klm_ii(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{
  
  if(multiplet_functions_ == 3) return;
  if (symmetry_) return;
  
  unsigned int i,j;
  
  unsigned int row_dim  =  ijk_.size();
  unsigned int col_dim =  begin_ij_;

  if (!row_dim || !col_dim) return;
  row_dim--; col_dim--;

  
#pragma omp parallel private(i,j)  
  {
    Blas_matrix block;
    
    for (i = 0; i <= row_dim + col_dim; i++ ) {
      
      unsigned int row_ind; 
      unsigned int adiag_size;
      unsigned int lower_dim = row_dim < col_dim ? row_dim : col_dim;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = row_dim + col_dim - i < lower_dim ?
	  row_dim + col_dim - i + 1  : lower_dim + 1; 
      }else {
	row_ind = i;
	adiag_size = i < lower_dim ? i + 1 : lower_dim + 1; 
      }
      unsigned int col_ind  = i - row_ind;
      
#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     
	
	unsigned int row = ijk_[row_ind - j];
	unsigned int col = col_ind + j;
      
	if (blocks_->block_klm_ii(configs_[row], configs_[col], block)) { 
	  
	  unsigned int block_rows = block.rows();
	  
	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0 , 1)
	      .daxpy((**vec_in)(col), block);
	    
	    (**vec_out)(col) +=  
	      block
	      * (**vec_in)(row, block_rows, 0, 1);
	    vec_out++; vec_in++;
	  }
	  
	}
      }    
    }
  }
}

// Builds the <3h1p,II|H|2h, open shell>
void Adc2_matrix::build_submatrix_klm_ij(Blas_matrix &adc2_matrix)
{
  Blas_matrix block;
  
  for(unsigned int row_inc, row = begin_ijk_; 
      row < configs_.size(); row += row_inc) {
    
    for(unsigned int col = begin_ij_; col < begin_jii_; col++) 
      
      if (blocks_->block_klm_ij(configs_[row], configs_[col], block))
	
	adc2_matrix(row, block.rows(), col , 1) = block;
    
    row_inc = multiplet_functions_
      *  blocks_->size_vir_group(configs_[row].vir);
    
  }
}

// Multiplies the <3h1p,II|H|2h, open shell> with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_klm_ij(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{
  unsigned int i,j;

  unsigned int row_dim  =  ijk_.size();
  unsigned int col_dim = begin_jii_ - begin_ij_;
  if (!row_dim || !col_dim) return;
  row_dim--; col_dim--;

#pragma omp parallel private(i,j)  
  {
    Blas_matrix block;
    for (i = 0; i <= row_dim + col_dim; i++ ) {
      
      unsigned int row_ind; 
      unsigned int adiag_size;
      unsigned int lower_dim = row_dim < col_dim ? row_dim : col_dim;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = row_dim + col_dim - i < lower_dim ? 
	  row_dim + col_dim - i + 1  : lower_dim + 1; 
      }else {
	row_ind = i;
	adiag_size = i < lower_dim ? i + 1 : lower_dim + 1; 
      }
      unsigned int col_ind  = i - row_ind;

#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     
	
	unsigned int row = ijk_[row_ind - j];
	unsigned int col =  begin_ij_ + col_ind + j;
      
	if (blocks_->block_klm_ij(configs_[row], configs_[col], block)) { 
	
	  unsigned int block_rows = block.rows();
	  
	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0 , 1)
	      .daxpy((**vec_in)(col), block);
	    
	    (**vec_out)(col) += 
	      block
	      * (**vec_in)(row, block_rows, 0, 1);
	    vec_out++; vec_in++;
	  }
	}
      }
    }
  }
}


// Builds the <3h1p,I|H|3h1p,I> block
void Adc2_matrix::build_submatrix_jii_lkk(Blas_matrix &adc2_matrix)
{
  Blas_matrix block;
  
  for (unsigned int row_inc, row = begin_jii_; 
       row < begin_ijk_; row += row_inc) {    
    
    for (unsigned int col_inc, col = begin_jii_; 
	 col <= row; col += col_inc) {
      
      if (blocks_->block_jii_lkk(configs_[row], configs_[col], block)) {

	adc2_matrix(row, block.rows(), col, block.cols()) = block;
      }
      
      col_inc =  blocks_->size_vir_group(configs_[col].vir);
    }
    
    row_inc =  blocks_->size_vir_group(configs_[row].vir);
  }
}



// Multiplies the <3h1p,I|H|3h1p,I> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_jii_lkk(Blas_matrix** block_out, 
					   Blas_matrix** block_in, int count)
{
  unsigned int i,j;
  // It's a square submatrix 
  unsigned int row_dim =  jii_.size();
  if (!row_dim) return;
  row_dim--; 

#pragma omp parallel private(i,j)  
  {
    
    Blas_matrix block;
    for (i = 0; i <= 2 * row_dim; i++ ) {

      unsigned int row_ind; 
      unsigned int adiag_size;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = (2 * row_dim - i) / 2 + 1;
      }else {
	row_ind = i;
	adiag_size = i / 2 + 1;
      }
      unsigned int col_ind  = i - row_ind;
      
#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     
	
	unsigned int row = jii_[row_ind - j];
	unsigned int col = jii_[col_ind + j];
	
	if (blocks_->block_jii_lkk(configs_[row], configs_[col], block)) { 
	  
	  unsigned int block_cols =  block.cols();		
	  unsigned int block_rows =  block.rows();

	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0, 1)
	      .dgemv('N', block, (**vec_in)(col, block_cols, 0, 1));
	    if (row != col) 
	      (**vec_out)(col, block_cols, 0, 1)
		.dgemv('T', block, (**vec_in)(row, block_rows, 0, 1));
	    vec_out++; vec_in++;
	  }
	}
      }
    }
  }
}

// Builds the <3h1p,II|H|3h1p,I> block
void Adc2_matrix::build_submatrix_ijk_mll(Blas_matrix &adc2_matrix)
{

  Blas_matrix block;

  for (unsigned int row_inc, row = begin_ijk_; 
       row < configs_.size(); row += row_inc) {    
    
    for (unsigned int col_inc, col = begin_jii_; 
	 col < begin_ijk_; col += col_inc) {
      
      if (blocks_->block_ijk_mll(configs_[row], configs_[col], block)) 
	
	adc2_matrix(row, block.rows(), col, block.cols()) = block;
	
      col_inc =  blocks_->size_vir_group(configs_[col].vir);
    }
    row_inc = multiplet_functions_
      *  blocks_->size_vir_group(configs_[row].vir);
  }
  
}

// Multiplies the <3h1p,II|H|3h1p,I> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_ijk_mll(Blas_matrix** block_out,
					 Blas_matrix** block_in, int count)
{

  unsigned int i,j;

  unsigned int row_dim  =  ijk_.size();
  unsigned int col_dim = jii_.size();
  if (!row_dim || !col_dim) return;
  row_dim--; col_dim--;


#pragma omp parallel private(i,j)  
  {
    Blas_matrix block;
    for (i = 0; i <= row_dim + col_dim; i++ ) {
      
      unsigned int row_ind; 
      unsigned int adiag_size;
      
      unsigned int lower_dim = row_dim < col_dim ? row_dim : col_dim;

      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = row_dim + col_dim - i < lower_dim ?
	  row_dim + col_dim - i + 1  : lower_dim + 1; 
      }else {
	row_ind = i;
	adiag_size = i < lower_dim ? i + 1 : lower_dim + 1; 
      }
      unsigned int col_ind  = i - row_ind;

#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     

	unsigned int row = ijk_[row_ind - j];
	unsigned int col = jii_[col_ind + j];
      
	if (blocks_->block_ijk_mll(configs_[row], configs_[col], block)) { 
	
	  unsigned int block_cols =  block.cols();		
	  unsigned int block_rows = block.rows();

	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0, 1)
	      .dgemv('N', block, (**vec_in)(col, block_cols, 0, 1));
	    (**vec_out)(col, block_cols, 0, 1)
	      .dgemv('T', block, (**vec_in)(row, block_rows, 0, 1));
	    vec_out++; vec_in++;
	  }
	}
      }
    }
  } //end omp parallel
}

// Builds the <3h1p,II|H|3h1p,II> block
void Adc2_matrix::build_submatrix_ijk_lmn(Blas_matrix &adc2_matrix)
{

  Blas_matrix block;
  
  for (unsigned int row_inc, row = begin_ijk_; 
       row < configs_.size(); row += row_inc) {    
    
    for (unsigned int col_inc, col = begin_ijk_; 
	 col <= row; col += col_inc) {
      
      if (blocks_->block_ijk_lmn(configs_[row], configs_[col], block)) 
	
	adc2_matrix(row, block.rows(), col, block.cols()) = block;
      
      col_inc = multiplet_functions_
	*  blocks_->size_vir_group(configs_[col].vir);
    }
    row_inc = multiplet_functions_
      *  blocks_->size_vir_group(configs_[row].vir);   
  }
}

// Multiplies the <3h1p,II|H|3h1p,II> block with the relevant part of the input vectors.
inline void Adc2_matrix::multiply_submatrix_ijk_lmn(Blas_matrix** block_out, 
					 Blas_matrix** block_in, int count)
{

  unsigned int i, j;
  unsigned int row_dim =  ijk_.size();
  if (!row_dim) return;
  row_dim--;

  
#pragma omp parallel private(i,j)  
  {
    Blas_matrix block;
    for (i = 0; i <= 2 * row_dim; i++ ) {

      unsigned int row_ind; 
      unsigned int adiag_size;
      
      if (i > row_dim) {
	row_ind = row_dim;
	adiag_size = (2 * row_dim - i) / 2 + 1;
      }else {
	row_ind = i;
	adiag_size = i / 2 + 1;
      }
      unsigned int col_ind  = i - row_ind;
      
#pragma omp for schedule(static, 1)
      for(j = 0; j < adiag_size; j++) {     

	unsigned int row = ijk_[row_ind - j];
	unsigned int col = ijk_[col_ind + j];
      
      
	if (blocks_->block_ijk_lmn(configs_[row], configs_[col], block)) { 
	
	  unsigned int block_cols = block.cols();
	  unsigned int block_rows = block.rows();
	  
	  Blas_matrix **vec_in = block_in; 
	  Blas_matrix **vec_out = block_out;
	  for(int vec = 0; vec < count; vec++) {
	    (**vec_out)(row, block_rows, 0, 1)
	      .dgemv('N', block, (**vec_in)(col, block_cols, 0, 1));
	    
	    if (row != col)
	      (**vec_out)(col, block_cols, 0, 1)
		.dgemv('T', block, (**vec_in)(row, block_rows, 0, 1));
	    
	    vec_out++; vec_in++;
	  }
	  
	}
      }
    }
  }
}



// Returns the configuration (a string) corresponding to the given row/column of 
// the double ionization adc2 matrix
string Adc2_matrix::get_conf(unsigned int i) const
{
  
    
  ostringstream o;
  if (i < begin_jii_)
    o << "<" << (int) configs_[i].occ[0]+1
      << "," << (int) configs_[i].occ[1]+1
      << "|";
  else if (i < begin_ijk_)
    o << "<" << (int) configs_[i].occ[0]+1
      << "," << (int) configs_[i].occ[1]+1
      << "," << (int) configs_[i].occ[2]+1
      << "," << configs_[i].vir+1
      << "|";
  else if (i < size()) {
    
    o << "<" << (int) configs_[i].occ[0]+1
      << "," << (int) configs_[i].occ[1]+1
      << "," << (int) configs_[i].occ[2]+1
      << "," << configs_[i].vir+1;
    switch (configs_[i].type) {
    case 0: 
      o << ",I"; break;
    case 1:
      o << ",II"; break;
    case 2: 
      o << ",III"; break;
    }
    o << "|";
  }  else
    o << "";

  return o.str();


    
}



void Adc2_matrix::print_configs()
{
  // print all configurations
  for(int i = 0; i < configs_.size(); i++) 
    cout << i << ' ' << get_conf(i) << endl;
    
}

#include "adc2_cap_matrix.hpp"
#include "zeroth_order.hpp"
#include "first_order.hpp"
#include "second_order.hpp"
#include "zeroth_order_triplet.hpp"
#include "first_order_triplet.hpp"
#include "second_order_triplet.hpp"
#include "../analysis/full_cap_analyzer.hpp"
#include "blas_matrix.hpp"


#include <iostream>
#include <sstream>
#include <string>
using namespace std;


#define SEC_
#define FIR_
//#undef SEC_
//#undef FIR_

Adc2_CAP_matrix::Adc2_CAP_matrix(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat,
				 unsigned int sym, unsigned int spin) :
  phis_(phis), symmetry_(sym), begin_ij_(0), begin_jii_(0), begin_ijk_(0), count_mult(0)
{

  switch (spin) {
    //Singlet
  case 0:
    multiplet_functions_ = 2;
    zero_ = new Zeroth_order_cap(phis_, tab, mat);
    fir_ = new First_order_cap(phis_, tab, mat);    
    sec_  = new Second_order_cap(phis_, tab, mat);
    break;
    //Triplet
  case 2:
    multiplet_functions_ = 3;
    zero_ = new Zeroth_order_cap_triplet(phis_, tab, mat);
    fir_ = new First_order_cap_triplet(phis_, tab, mat);    
    sec_  = new Second_order_cap_triplet(phis_, tab, mat);
    break;
  default:
    throw string("Only singlet ADC2 double ionization CAP is implemented for now...\n");
  }
  
  // Determine the configurations for the given symmetry and spin.
  add_ii_configs();
  add_ij_configs();
  add_jiir_configs();
  add_ijkr_configs();
  
  
    print_configs();
  
}

int Adc2_CAP_matrix::accept_analyzer(ADC_analyzer& analyst)
{
  return analyst.analyze(*this);
}

Adc2_CAP_matrix::~Adc2_CAP_matrix()
{
  delete zero_;
  delete sec_;
  delete fir_;
}


// Adds closed shell 2h configurations, |ii>
void Adc2_CAP_matrix::add_ii_configs()
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
void Adc2_CAP_matrix::add_ij_configs()
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
void Adc2_CAP_matrix::add_jiir_configs()
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
void Adc2_CAP_matrix::add_ijkr_configs()
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
void Adc2_CAP_matrix::build_matrix(Blas_matrix &adc2_matrix)
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

int Adc2_CAP_matrix::operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count)
{  
  
  if (count_mult < 2) {

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

  
  count_mult ++;

  return count;

}


// Builds the <2h, closed shell|H|2h, closed shell> block
void Adc2_CAP_matrix::build_submatrix_ii_jj(Blas_matrix &adc2_matrix)
{
  
  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;
  

  
  for(unsigned int row = 0; row < begin_ij_; row++) 
    for(unsigned int col = 0; col < begin_ij_; col++) {
      double element = 0.;
      zero_->block_ii_jj(configs_[row], configs_[col], element);
#ifdef SEC_
      sec_->block_ii_jj(configs_[row], configs_[col], element);
#endif
      adc2_matrix(row, col) = element;
      
    }

}

// Multiplies the <2h, closed shell|H|2h, closed shell> block with the relevant part of the input vectors.
inline void 
Adc2_CAP_matrix::multiply_submatrix_ii_jj(Blas_matrix** block_out, 
					  Blas_matrix** block_in, int count)
{
  

}


//Builds the <2h, open shell|H|2h, closed shell> block
void Adc2_CAP_matrix::build_submatrix_ij_kk(Blas_matrix &adc2_matrix)
{
  
  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;
  


  for(unsigned int row = begin_ij_; row < begin_jii_; row++) 
    for(unsigned int col = 0; col < begin_ij_; col++) {
      double element = 0.;
      zero_->block_ij_kk(configs_[row], configs_[col], element);
#ifdef SEC_
      sec_->block_ij_kk(configs_[row], configs_[col], element);
#endif
      adc2_matrix(row, col) = adc2_matrix(col, row) = element;
      
    }
}

// Multiplies the <2h, open shell|H|2h, closed shell> block  with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_ij_kk(Blas_matrix** block_out, 
				       Blas_matrix** block_in, int count)
{
}

// Builds the <2h, open shell|H|2h, open shell> block
void Adc2_CAP_matrix::build_submatrix_ij_kl(Blas_matrix &adc2_matrix)
{

  
  for(unsigned int row = begin_ij_; row < begin_jii_; row++) 
    for(unsigned int col = begin_ij_; col < begin_jii_; col++) {
      double element = 0.;
      zero_->block_ij_kl(configs_[row], configs_[col], element);
#ifdef SEC_
      sec_->block_ij_kl(configs_[row], configs_[col], element);
#endif
      adc2_matrix(row, col) = element;
      
    }
}


// Multiplies the <2h, open shell|H|2h, open shell> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_ij_kl(Blas_matrix** block_out, 
					 Blas_matrix** block_in, int count)
{

}

// Builds the <3h1p,I|H|2h, closed shell> block
void Adc2_CAP_matrix::build_submatrix_lkk_ii(Blas_matrix &adc2_matrix)
{

  if(multiplet_functions_ == 3) return;
  if (symmetry_) return;

  for(unsigned int row_inc, row = begin_jii_; 
      row < begin_ijk_; row += row_inc) {
    
    for(unsigned int col = 0; col < begin_ij_; col++) {
      
      Blas_matrix block;
      bool z = zero_->block_lkk_ii(configs_[row], configs_[col], block);
#ifndef FIR_
      bool f = false;
#endif
#ifdef FIR_
      bool f = fir_->block_lkk_ii(configs_[row], configs_[col], block);
#endif
      if (z || f) {
	//block *= -1.;
	adc2_matrix(row, block.rows(), col, 1) =  block;
      }
    }
    row_inc =  zero_->size_vir_group(configs_[row].vir);   
  }  
}


// Multiplies the <3h1p,I|H|2h, closed shell> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_lkk_ii(Blas_matrix** block_out, 
					  Blas_matrix** block_in, int count)
{

}


// Builds the <3h1p,I|H|2h, open shell> block
void Adc2_CAP_matrix::build_submatrix_lkk_ij(Blas_matrix &adc2_matrix)
{

  for(unsigned int row_inc, row = begin_jii_; 
      row < begin_ijk_; row += row_inc) {
    
    for(unsigned int col = begin_ij_; col < begin_jii_; col++) {
      Blas_matrix block;
      bool z = zero_->block_lkk_ij(configs_[row], configs_[col], block);
#ifndef FIR_
      bool f = false;
#endif
#ifdef FIR_
      bool f = fir_->block_lkk_ij(configs_[row], configs_[col], block);
#endif
      if (z || f) {
	//block *= -1.;
	adc2_matrix(row, block.rows(), col, 1) =  block;
      }
    }
      
    
    row_inc =  zero_->size_vir_group(configs_[row].vir);
  }  
}


// Multiplies the <3h1p,I|H|2h, open shell> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_lkk_ij(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{

}

// Builds the <3h1p,II|H|2h, closed shell>
void Adc2_CAP_matrix::build_submatrix_klm_ii(Blas_matrix &adc2_matrix)
{

  if(multiplet_functions_ == 3) return;
  if(symmetry_) return;

  for(unsigned int row_inc, row = begin_ijk_; 
      row < configs_.size(); row += row_inc) {
    
    for(unsigned int col = 0; col < begin_ij_; col++) {
      Blas_matrix block;
      bool z = zero_->block_klm_ii(configs_[row], configs_[col], block);
#ifndef FIR_
      bool f = false;
#endif
#ifdef FIR_
      bool f = fir_->block_klm_ii(configs_[row], configs_[col], block);
#endif
      if (z || f) {
	//block *= -1.;
	adc2_matrix(row, block.rows(), col, 1) =  block;
      }
    }

    row_inc = multiplet_functions_
      *  zero_->size_vir_group(configs_[row].vir);
    
  }
}

// Multiplies the <3h1p,II|H|2h, closed shell> with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_klm_ii(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{
  

}

// Builds the <3h1p,II|H|2h, open shell>
void Adc2_CAP_matrix::build_submatrix_klm_ij(Blas_matrix &adc2_matrix)
{
  
  for(unsigned int row_inc, row = begin_ijk_; 
      row < configs_.size(); row += row_inc) {
    
    for(unsigned int col = begin_ij_; col < begin_jii_; col++)  {
      Blas_matrix block;
      
      bool z = zero_->block_klm_ij(configs_[row], configs_[col], block);
#ifndef FIR_
      bool f = false;
#endif
#ifdef FIR_
      bool f = fir_->block_klm_ij(configs_[row], configs_[col], block);
#endif
      if (z || f) {
	//block *= -1.;
	adc2_matrix(row, block.rows(), col, 1) =  block;
      }
    }
    
    row_inc = multiplet_functions_
      *  zero_->size_vir_group(configs_[row].vir);
    
  }
}

// Multiplies the <3h1p,II|H|2h, open shell> with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_klm_ij(Blas_matrix** block_out, 
					Blas_matrix** block_in, int count)
{

}


// Builds the <3h1p,I|H|3h1p,I> block
void Adc2_CAP_matrix::build_submatrix_jii_lkk(Blas_matrix &adc2_matrix)
{
  Blas_matrix block;
  
  for (unsigned int row_inc, row = begin_jii_; 
       row < begin_ijk_; row += row_inc) {    
    
    for (unsigned int col_inc, col = begin_jii_; 
	 col <= row; col += col_inc) {
      
      if (zero_->block_jii_lkk(configs_[row], configs_[col], block)) {

	adc2_matrix(row, block.rows(), col, block.cols()) = block;
      }
      
      col_inc =  zero_->size_vir_group(configs_[col].vir);
    }
    
    row_inc =  zero_->size_vir_group(configs_[row].vir);
  }
}



// Multiplies the <3h1p,I|H|3h1p,I> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_jii_lkk(Blas_matrix** block_out, 
					   Blas_matrix** block_in, int count)
{

  
}
void Adc2_CAP_matrix::build_submatrix_ijk_mll(Blas_matrix &adc2_matrix)
{

  Blas_matrix block;

  for (unsigned int row_inc, row = begin_ijk_;
       row < configs_.size(); row += row_inc) {

    for (unsigned int col_inc, col = begin_jii_;
         col < begin_ijk_; col += col_inc) {

      if (zero_->block_ijk_mll(configs_[row], configs_[col], block)) {
	//block *= -1.;
        adc2_matrix(row, block.rows(), col, block.cols()) = block;
      }

      col_inc =  zero_->size_vir_group(configs_[col].vir);
    }
    row_inc = multiplet_functions_
      *  zero_->size_vir_group(configs_[row].vir);
  }

}

// Multiplies the <3h1p,II|H|3h1p,I> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_ijk_mll(Blas_matrix** block_out,
					 Blas_matrix** block_in, int count)
{

}

// Builds the <3h1p,II|H|3h1p,II> block
void Adc2_CAP_matrix::build_submatrix_ijk_lmn(Blas_matrix &adc2_matrix)
{

  Blas_matrix block;
  
  for (unsigned int row_inc, row = begin_ijk_; 
       row < configs_.size(); row += row_inc) {    
    
    for (unsigned int col_inc, col = begin_ijk_; 
	 col <= row; col += col_inc) {
      
      if (zero_->block_ijk_lmn(configs_[row], configs_[col], block)) 
	
	adc2_matrix(row, block.rows(), col, block.cols()) = block;
      
      col_inc = multiplet_functions_
	*  zero_->size_vir_group(configs_[col].vir);
    }
    row_inc = multiplet_functions_
      *  zero_->size_vir_group(configs_[row].vir);   
  }
}

// Multiplies the <3h1p,II|H|3h1p,II> block with the relevant part of the input vectors.
inline void Adc2_CAP_matrix::multiply_submatrix_ijk_lmn(Blas_matrix** block_out, 
					 Blas_matrix** block_in, int count)
{

}



// Returns the configuration (a string) corresponding to the given row/column of 
// the double ionization adc2 matrix
string Adc2_CAP_matrix::get_conf(unsigned int i) const
{
  
    
  ostringstream o;
  if (i < begin_jii_)
    o << '<' << (int) configs_[i].occ[0]+1
      << ',' << (int) configs_[i].occ[1]+1
      << '|';
  else if (i < begin_ijk_)
    o << '<' << (int) configs_[i].occ[0]+1
      << ',' << (int) configs_[i].occ[1]+1
      << ',' << (int) configs_[i].occ[2]+1
      << ',' << configs_[i].vir+1
      << '|';
  else if (i < size()) {
    
    o << '<' << (int) configs_[i].occ[0]+1
      << ',' << (int) configs_[i].occ[1]+1
      << ',' << (int) configs_[i].occ[2]+1
      << ',' << configs_[i].vir+1;
    switch (configs_[i].type) {
    case 0: 
      o << ",I"; break;
    case 1:
      o << ",II"; break;
    case 2: 
      o << ",III"; break;
    }
    o << '|';
  }  else
    o << "";

  return o.str();


    
}



void Adc2_CAP_matrix::print_configs()
{
  // print all configurations
  for(int i = 0; i < configs_.size(); i++) 
    cout << i << ' ' << get_conf(i) << endl;
    
}

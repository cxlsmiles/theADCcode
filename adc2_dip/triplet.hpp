#ifndef __TRIPLET_HPP__
#define __TRIPLET_HPP__


// The file contains the declaration of the Triplet class which implements 
// the triplet formulas for the building blocks of the ADC2 matrix.
// See Chem. Phis. 329, 11: Tables A.2, A.3, A.5, A.7
// and eqs. A.1, A.2, A.3. There are misprints
// in the formulas!

#include "adc2_dip_blocks.hpp"

class Integral_table;
class Triplet: public ADC2_DIP_blocks {
  
  double w_term(unsigned int i, unsigned int j,
		unsigned int s, unsigned int r);
  double u_term(unsigned int i, unsigned int j, 
		unsigned int k, unsigned int l,
		unsigned int s, unsigned int r);

public:


  Triplet(SCF_data_reader& phis, Integral_table& tab) : ADC2_DIP_blocks(phis,tab) {}


  // Closed shell blocks do not exist for triplet
  bool block_ii_jj(const Config& row, const Config& col, double &element)     {return false;}
  bool block_ij_kk(const Config& row, const Config& col, double &element)     {return false;}
  bool block_lkk_ii(const Config& row, const Config& col, Blas_matrix& block) {return false;}
  bool block_klm_ii(const Config& row, const Config& col, Blas_matrix& block) {return false;}

  // Open shell blocks.
  bool block_ij_kl(const Config& row, const Config& col, double &element);
  bool block_lkk_ij(const Config& row, const Config& col, Blas_matrix& block);
  bool block_klm_ij(const Config& row, const Config& col, Blas_matrix& block);
  bool block_jii_lkk(const Config& row, const Config& col, Blas_matrix& block);
  bool block_ijk_mll(const Config& row, const Config& col, Blas_matrix& block);
  bool block_ijk_lmn(const Config& row, const Config& col, Blas_matrix& block);

};

#endif //#ifndef __TRIPLET_HPP__

#ifndef __SINGLET_HPP__
#define __SINGLET_HPP__

// The file contains the declaration of the Singlet class which implements 
// the singlet formulas for the building blocks of the ADC2 matrix.
// See Chem. Phis. 329, 11: Tables A.2, A.3, A.4, A.6
// and eqs. A.1, A.2, A.3, A.4, A.5. There are misprints
// in the formulas!

#include "adc2_dip_blocks.hpp"

class Integral_table;
class Singlet: public ADC2_DIP_blocks {
  
  double w_term(unsigned i, unsigned j, unsigned int s, unsigned int r);
  double u_term(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int s, unsigned int r);

public:


  Singlet(SCF_data_reader& phis, Integral_table& tab) : ADC2_DIP_blocks(phis,tab) {}
  
  bool block_ii_jj(const Config& row, const Config& col, double &element);
  bool block_ij_kk(const Config& row, const Config& col, double &element);
  bool block_ij_kl(const Config& row, const Config& col, double &element);
  bool block_lkk_ii(const Config& row, const Config& col, Blas_matrix& block);
  bool block_lkk_ij(const Config& row, const Config& col, Blas_matrix& block);
  bool block_klm_ii(const Config& row, const Config& col, Blas_matrix& block);
  bool block_klm_ij(const Config& row, const Config& col, Blas_matrix& block);
  bool block_jii_lkk(const Config& row, const Config& col, Blas_matrix& block);
  bool block_ijk_mll(const Config& row, const Config& col, Blas_matrix& block);
  bool block_ijk_lmn(const Config& row, const Config& col, Blas_matrix& block);
  
};

#endif //#ifndef __SINGLET_HPP__

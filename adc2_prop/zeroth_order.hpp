#ifndef __ZEROTH_ORDER_HPP__
#define __ZEROTH_ORDER_HPP__

#include "scf_data/scf_data_reader.hpp"
#include "blas_matrix.hpp"
#include <vector>
#include "matrices.hpp"
#include "adc2_dip/adc2_dip_blocks.hpp"

struct Config;
class Integral_table;
class Zeroth_order_cap:  public ADC2_DIP_blocks {

  SCF_data_reader& phis_;
  Blas_matrix Doo;
  std::vector<Blas_matrix> Dvo;
  std::vector<Blas_matrix> Dvv;
  std::vector<unsigned int> vir_group_sizes_;  
  // Stores the number of virtual orbitals in a given symmetry.

  double gs_expt_value_;
  
  inline unsigned sym(unsigned orb)  
  { return phis_.irrep(orb); } // get the symmetry of an orbital


public:


  Zeroth_order_cap(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat);

  inline unsigned size_vir_group(unsigned vir) 
  {return vir_group_sizes_[phis_.irrep(vir)];}  
  // Get the number of virtual orbitals for a given symmetry

  
  virtual bool block_ii_jj(const Config& row, const Config& col, double &element);
  virtual bool block_ij_kk(const Config& row, const Config& col, double &element);
  virtual bool block_ij_kl(const Config& row, const Config& col, double &element);
  virtual bool block_lkk_ii(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_lkk_ij(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_klm_ii(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_klm_ij(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_jii_lkk(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_ijk_mll(const Config& row, const Config& col, Blas_matrix& block);
  virtual bool block_ijk_lmn(const Config& row, const Config& col, Blas_matrix& block);
  
};

#endif //#ifndef __ZEROTH_ORDER_HPP__

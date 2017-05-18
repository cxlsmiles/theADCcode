#ifndef __SECOND_ORDER_HPP__
#define __SECOND_ORDER_HPP__

#include "scf_data/scf_data_reader.hpp"
#include "adc2_dip/adc2_dip_blocks.hpp"
#include <vector>
#include "matrices.hpp"

struct Config;
class Blas_matrix;
class Integral_table;

class Second_order_cap:  public ADC2_DIP_blocks {

  SCF_data_reader& phis_;

  Triangular_matrix<double> d;

  Blas_matrix* rho;

  std::vector<unsigned> occs_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> virs_[8];
  double gs_expt_value2_;


  double V1212(unsigned a, unsigned b, unsigned c, unsigned d);
  double d_ij_i1j1(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_1(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_2(unsigned i, unsigned j, unsigned i1, unsigned j1);  
  double d2_3(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_4(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_5(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_6(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_7(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_8(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_9(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_10(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_11(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_12(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_13(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_14(unsigned i, unsigned j, unsigned i1, unsigned j1);
  double d2_15(unsigned i, unsigned j, unsigned i1, unsigned j1);

  double vir_sum_rho_d(unsigned i, unsigned i1);


public:


  Second_order_cap(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat);

  
  virtual bool block_ii_jj(const Config& row, const Config& col, double &element);
  virtual bool block_ij_kk(const Config& row, const Config& col, double &element);
  virtual bool block_ij_kl(const Config& row, const Config& col, double &element);
  

  
  virtual bool block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }

  virtual bool block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }
  
  virtual bool block_klm_ii(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }

  virtual bool block_klm_ij(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }

  virtual bool block_jii_lkk(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }


  virtual bool block_ijk_mll(const Config &row, const Config &col, Blas_matrix& block)
  {
    return false;
  }
  

  virtual bool block_ijk_lmn(const Config &row, const Config &col, Blas_matrix &block)
  {
    return false;
  }



};

#endif //#ifndef __SECOND_ORDER_HPP__

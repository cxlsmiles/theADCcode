#ifndef __FIRST_ORDER_HPP__
#define __FIRST_ORDER_HPP__

#include "scf_data/scf_data_reader.hpp"
#include "adc2_dip/adc2_dip_blocks.hpp"
#include <vector>
#include "matrices.hpp"

struct Config;
class Integral_table;

class First_order_cap:  public ADC2_DIP_blocks {

  SCF_data_reader& phis_;



  double V1212(unsigned a, unsigned b, unsigned c, unsigned d);
  Triangular_matrix<double> d;


  std::vector<unsigned int> vir_group_sizes_;  
  double bk_sum(unsigned a1, unsigned x1);
  double b_sum(unsigned a1, unsigned x1, unsigned y1, unsigned z1);

  std::vector<unsigned> occs_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> virs_[8];


  inline unsigned size_vir_group(unsigned vir) 
  {return vir_group_sizes_[phis_.irrep(vir)];}  
  // Get the number of virtual orbitals for a given symmetry

  inline unsigned sym(unsigned orb)  
  { return phis_.irrep(orb); } // get the symmetry of an orbital


public:


  First_order_cap(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat);

  
  bool block_ii_jj(const Config& row, const Config& col, double &element)
  { return false;}
  bool block_ij_kk(const Config& row, const Config& col, double &element)
  { return false;}
  bool block_ij_kl(const Config& row, const Config& col, double &element)
  { return false;}
  


  virtual bool block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block);


  virtual bool block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block);


  virtual bool block_klm_ii(const Config &row, const Config &col, Blas_matrix& block);


  virtual bool block_klm_ij(const Config &row, const Config &col, Blas_matrix& block);



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

#endif //#ifndef __FIRST_ORDER_HPP__

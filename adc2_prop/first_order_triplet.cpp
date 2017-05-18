#include "first_order_triplet.hpp"
#include "adc2_dip/config.hpp"
#include "matrices.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

inline double First_order_cap_triplet::V1212(unsigned a, unsigned b, unsigned c, unsigned d)
{
  return V1122(a,c,b,d);
}


First_order_cap_triplet::First_order_cap_triplet(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat) 
  : phis_(phis), ADC2_DIP_blocks(phis, tab), d(mat->rows())
{


  for(unsigned row = 0; row < mat->rows(); row++)
    for(unsigned col = 0; col <= row; col++) {
      d(row, col) = (*mat)(row,col);
    }
  
  //Group the occupied and virtual orbitals according to symmetry;
  for(unsigned i = 0; i < phis_.number_occupied();  i++) 
    occs_[phis_.irrep(i)].push_back(i); 

  for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals();  i++)
    virs_[phis_.irrep(i)].push_back(i); 



  std::vector<unsigned int> sym_block_location_;
  
  vir_group_sizes_.assign(phis_.number_irreps(), 0);
  sym_block_location_.assign(phis_.number_orbitals(),0);
  
  
  for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals(); i++) 
    sym_block_location_[i] = vir_group_sizes_[phis_.irrep(i)]++;
  


}


inline double First_order_cap_triplet::b_sum(unsigned a1, unsigned x1, unsigned y1, unsigned z1)
{
  
  unsigned a1_sym = phis_.irrep(a1);
  unsigned x1_sym = phis_.irrep(x1);
  unsigned y1_sym = phis_.irrep(y1);
  double sum_b = 0.;
  unsigned b_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, x1_sym), y1_sym);
  if (b_sym != phis_.irrep(z1)) return 0.;
  for(unsigned b_  = 0; b_ < virs_[b_sym].size(); b_++) {
    unsigned b = virs_[b_sym][b_];

    sum_b += V1212(a1,b,x1,y1) * d(b,z1)
      / (phis_.energy(a1)+phis_.energy(b)-phis_.energy(x1)-phis_.energy(y1));
    
  }
  

  return sum_b;
}


inline double First_order_cap_triplet::bk_sum(unsigned a1, unsigned x1)
{
  
  unsigned a1_sym = phis_.irrep(a1);
  unsigned x1_sym = phis_.irrep(x1);
  double sum_bk = 0.;
  for(unsigned b  = phis_.number_occupied(); b < phis_.number_orbitals(); b++) {
    unsigned k_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, phis_.irrep(b)), x1_sym);
    for(unsigned k_  = 0; k_ < occs_[k_sym].size(); k_++) {
      unsigned k = occs_[k_sym][k_];
      sum_bk += (2. * V1212(a1,b,x1,k) - V1212(a1,b,k,x1)) * d(b,k)
	/ (phis_.energy(a1)+phis_.energy(b)-phis_.energy(x1)-phis_.energy(k));
      
    }
  }

  return sum_bk;
}


bool First_order_cap_triplet::block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  return false;
}

bool First_order_cap_triplet::block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned i1 = row.occ[0];
  unsigned j1 = row.occ[1];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_ii1 = i == i1;
  bool delta_ij1 = i == j1;
  bool delta_jj1 = j == j1;
  bool delta_ji1 = j == i1;

  if (!delta_ii1 && !delta_jj1 && !delta_ij1 && !delta_ji1)
    return false;


  unsigned int block_size = size_vir_group(row.vir);
  unsigned vir_sym = sym(row.vir);

  if (!block.size()) {
    block.allocate(block_size);
    block = 0.;
  }

  if (delta_ii1 && delta_jj1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += bk_sum(vir,j1);
    }
  }

  if (delta_ij1 && delta_ji1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) -= bk_sum(vir,j1);
    }
  }


  if (delta_ii1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += -b_sum(vir,j1,j1,j);
    }
  }

  if (delta_ji1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += b_sum(vir,j1,j1,i);
    }
  }
  if (delta_ij1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += b_sum(vir,j1,i1,j);
    }
  }
  if (delta_jj1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += -b_sum(vir,j1,i1,i);
    }
  }



  return true;
}


bool First_order_cap_triplet::block_klm_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  return false;
}


bool First_order_cap_triplet::block_klm_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned i1 = row.occ[0];
  unsigned j1 = row.occ[1];
  unsigned k1 = row.occ[2];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_ii1 = i == i1;
  bool delta_ij1 = i == j1;
  bool delta_ik1 = i == k1;
  bool delta_ji1 = j == i1;
  bool delta_jj1 = j == j1;
  bool delta_jk1 = j == k1;

  if (!delta_ii1 && !delta_ij1 && !delta_ik1 && !delta_ji1 && !delta_jj1 && !delta_jk1)
    return false;


  unsigned int block_size = size_vir_group(row.vir);
  unsigned vir_sym = sym(row.vir);
  
  if (!block.size()) {
    block.allocate(3 * block_size);
    block = 0.;
  }


  if (delta_ii1 && delta_jj1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];

      //block(vir_) += 0;

      double sum =  bk_sum(vir,k1);
      //block0
      
      block(vir_+block_size) -= sum;
      
      block(vir_+2*block_size) -= sum;

    }
  }

  if (delta_ij1 && delta_jk1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];

      double sum = bk_sum(vir,i1);

      block(vir_) -= sum;
      
      block(vir_+block_size) -= sum;
      //block2
    }
  }

  if (delta_ii1 && delta_jk1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      double sum = bk_sum(vir,j1);

      block(vir_) += sum;
      //block1
      block(vir_+2*block_size) += sum;

    }
  }  
  

  if (delta_ii1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += -b_sum(vir,j1,k1,j);
      block(vir_+block_size) += b_sum(vir,k1,j1,j);
      block(vir_+2*block_size) += b_sum(vir,k1,j1,j)-b_sum(vir,j1,k1,j);
    }
  }
  if (delta_ji1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += b_sum(vir,j1,k1,i);
      block(vir_+block_size) += -b_sum(vir,k1,j1,i);
      block(vir_+2*block_size) += -b_sum(vir,k1,j1,i)+b_sum(vir,j1,k1,i);
    }
  }
  if (delta_ij1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += b_sum(vir,i1,k1,j);
      block(vir_+block_size) += -b_sum(vir,k1,i1,j)+b_sum(vir,i1,k1,j);
      block(vir_+2*block_size) += -b_sum(vir,k1,i1,j);
    }
  }
  if (delta_jj1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += -b_sum(vir,i1,k1,i);
      block(vir_+block_size) += b_sum(vir,k1,i1,i)-b_sum(vir,i1,k1,i);
      block(vir_+2*block_size) += b_sum(vir,k1,i1,i);
    }
  }
  if (delta_ik1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += b_sum(vir,j1,i1,j)-b_sum(vir,i1,j1,j);
      block(vir_+block_size) += -b_sum(vir,i1,j1,j);
      block(vir_+2*block_size) += b_sum(vir,j1,i1,j);
    }
  }
  if (delta_jk1) {
    for(unsigned vir_ = 0; vir_ < virs_[vir_sym].size(); vir_++) {
      unsigned vir = virs_[vir_sym][vir_];
      block(vir_) += -b_sum(vir,j1,i1,i)+b_sum(vir,i1,j1,i);
      block(vir_+block_size) += b_sum(vir,i1,j1,i);
//      block(vir_+block_size) += -b_sum(vir,j1,i1,i);
// NICO 24.03.2011
      block(vir_+2*block_size) += -b_sum(vir,j1,i1,i);
    }
  }



  return true;
}


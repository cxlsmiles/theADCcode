#include "zeroth_order_triplet.hpp"
#include "adc2_dip/config.hpp"
#include "matrices.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

// JUST A TEST!

Zeroth_order_cap_triplet::Zeroth_order_cap_triplet(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat) 
  : phis_(phis), gs_expt_value_(0.), ADC2_DIP_blocks(phis, tab)
{
  
  Doo.allocate(phis_.number_occupied(), phis_.number_occupied()); 
  for(unsigned occ_row = 0; occ_row<phis_.number_occupied(); occ_row++) 
    for(unsigned occ_col = 0; occ_col<phis_.number_occupied(); occ_col++)
      //Doo(occ_row, occ_col) = 10 *  occ_row + occ_col;
      Doo(occ_row, occ_col) = (*mat)(occ_row,occ_col);

  for(unsigned occ_row = 0; occ_row<phis_.number_occupied(); occ_row++) 
   gs_expt_value_ += Doo(occ_row, occ_row);
  
  gs_expt_value_ *= 2;

  std::vector<unsigned int> sym_block_location_;
  
  vir_group_sizes_.assign(phis_.number_irreps(), 0);
  sym_block_location_.assign(phis_.number_orbitals(),0);
  
  
  for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals(); i++) 
    sym_block_location_[i] = vir_group_sizes_[phis_.irrep(i)]++;
  
  
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) {
    Dvv.push_back(Blas_matrix(vir_group_sizes_[sym], vir_group_sizes_[sym]));
    Dvo.push_back(Blas_matrix(vir_group_sizes_[sym], phis_.number_occupied()));
    
    for(unsigned int vir_row = phis_.number_occupied(); vir_row < phis_.number_orbitals(); vir_row++){
      if (phis_.irrep(vir_row) == sym)
	for(unsigned occ = 0; occ<phis_.number_occupied(); occ++)
	  //Dvo[sym](sym_block_location_[vir_row], occ) = 100 * vir_row + occ;
	  Dvo[sym](sym_block_location_[vir_row], occ) = (*mat)(vir_row,occ);
      for(unsigned int vir_col = phis_.number_occupied(); vir_col < phis_.number_orbitals(); vir_col++)
	if ((phis_.irrep(vir_row) == sym) && (phis_.irrep(vir_col) == sym))
// 	  Dvv[sym](sym_block_location_[vir_row], sym_block_location_[vir_col]) 
// 	    = 100 * vir_row + vir_col;
	  Dvv[sym](sym_block_location_[vir_row], sym_block_location_[vir_col]) 
 	    =  (*mat)(vir_row,vir_col);
    }
  }


}


bool Zeroth_order_cap_triplet::block_ii_jj(const Config &row, const Config &col, double &element)
{
  return false;
  
}


bool Zeroth_order_cap_triplet::block_ij_kk(const Config &row, const Config &col, double &element)
{

  return false;
  
}

bool Zeroth_order_cap_triplet::block_ij_kl(const Config &row, const Config &col, double &element)
{

  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  
  unsigned k = col.occ[0];
  unsigned l = col.occ[1];
  
  bool delta_ik = i == k;
  bool delta_jk = j == k;
  bool delta_il = i == l;
  bool delta_jl = j == l;

  
  element = 0.;

  if (delta_jl && delta_ik)
    element += gs_expt_value_;

  if (delta_jl)
    element -= Doo(k,i);
  if (delta_ik)
    element -= Doo(l,j);
  if (delta_jk)
    element += Doo(l,i);
  if (delta_il)
    element += Doo(k,j);

  return true;
  
}


bool Zeroth_order_cap_triplet::block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  return false;
}

bool Zeroth_order_cap_triplet::block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned l = row.occ[0];
  unsigned k = row.occ[1];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_jk = j == k;
  bool delta_il = i == l;
  
  bool delta_ik = i == k;
  bool delta_jl = j == l;


  if (!(delta_il && delta_jk) && !(delta_ik && delta_jl))
    return false;
  
  unsigned int block_size = size_vir_group(row.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(block_size);
  
  block = 0.;
  Submatrix submat = block(0,block_size, 0, 1);
  
  if (delta_il && delta_jk)
    submat.daxpy(1., Dvo[sym_vir](0, block_size, k, 1));  

  if (delta_ik && delta_jl)
    submat.daxpy(-1., Dvo[sym_vir](0, block_size, k, 1));  


  return true;
}


bool Zeroth_order_cap_triplet::block_klm_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  return false;
}


bool Zeroth_order_cap_triplet::block_klm_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned k = row.occ[0];
  unsigned l = row.occ[1];
  unsigned m = row.occ[2];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_ik = i == k;
  bool delta_il = i == l;
  bool delta_jl = j == l;
  bool delta_jm = j == m;
  
  if (!(delta_jm && delta_ik) && !(delta_jm && delta_il)
      && !(delta_jl && delta_ik))
    return false;
  
  unsigned int block_size = size_vir_group(row.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(3 * block_size);
  block = 0.;

  Submatrix part[] = {block(0, block_size, 0, 1),
		      block(block_size, block_size, 0, 1),
		      block(2*block_size, block_size, 0, 1)};

  if (delta_jm && delta_il) {
    
    part[0].daxpy(-1., Dvo[sym_vir](0, block_size, k, 1));
    part[1].daxpy(-1., Dvo[sym_vir](0, block_size, k, 1));
    //part2
  }

  if (delta_jm && delta_ik) {
    part[0].daxpy(1., Dvo[sym_vir](0, block_size, l, 1));
    //part1
    part[2].daxpy(1., Dvo[sym_vir](0, block_size, l, 1));
  }
  
  
  if (delta_jl && delta_ik) {
    //part0
    part[1].daxpy(-1., Dvo[sym_vir](0, block_size, m, 1));
    part[2].daxpy(-1., Dvo[sym_vir](0, block_size, m, 1));
  }


  return true;
}


bool Zeroth_order_cap_triplet::block_jii_lkk(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned j = row.occ[0];
  unsigned i = row.occ[1];
  
  unsigned l = col.occ[0];
  unsigned k = col.occ[1];
  
 
  bool delta_jl = j == l;
  bool delta_ik = i == k;

  bool delta_jk = j == k;
  bool delta_il = i == l;

  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  if (!delta_sym) return false;

  if (
      !(delta_jl && delta_ik) &&
      !(delta_ik)             &&
      !(delta_jk && delta_il) 
      )
    return false;

  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(row_block_size, col_block_size);
  block = 0.;

  double diag_term = 0.;  

  if (delta_jl && delta_ik) {

    block.daxpy(1., Dvv[sym_vir]);
    diag_term += gs_expt_value_ - 2. * Doo(k,i);
  }

  if (delta_ik) 
    diag_term -= Doo(l,j);


  if (delta_jk && delta_il) 
    diag_term += Doo(k,i);

  block.add_diag(diag_term);

  
  return true;
}


bool Zeroth_order_cap_triplet::block_ijk_mll(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  unsigned k = row.occ[2];
  
  unsigned m = col.occ[0];
  unsigned l = col.occ[1];
  

  bool delta_jm = j == m;
  bool delta_im = i == m;
  bool delta_kl = k == l;
  bool delta_jl = j == l;

  bool delta_il = i == l;
  bool delta_km = k == m;

  
  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  
  if (!delta_sym) return false;

  if (
      !(delta_jm && delta_kl) &&
      !(delta_im && delta_jl) &&
      !(delta_im && delta_kl) &&
      !(delta_il && delta_jm) &&
      !(delta_il && delta_km) &&
      !(delta_jl && delta_km) 
      )
    return false;

  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(3 * row_block_size, col_block_size);
  block = 0.;
  
  Submatrix part[] = 
    {block(0, row_block_size, 0, col_block_size),
     block(row_block_size, row_block_size, 0, col_block_size),
     block(2*row_block_size, row_block_size, 0, col_block_size)};

  double diag_term[] = {0., 0., 0.};
  
    
  if (delta_jm && delta_kl) {
    diag_term[0] += Doo(l,i);
    //term1
    diag_term[2] -= Doo(l,i);
  }
  if (delta_im && delta_jl) {
    diag_term[0] -= Doo(l,k);
    diag_term[1] += Doo(l,k);
    //term2
  }

  if (delta_im && delta_kl) {
    diag_term[0] -= Doo(l,j);
    diag_term[1] += Doo(l,j);
    //term2
  }

  if (delta_il && delta_jm) {
    diag_term[0] += Doo(l,k);
    //term1
    diag_term[2] -= Doo(l,k);
  }
  if (delta_il && delta_km) {
    //term0 
    diag_term[1] -= Doo(l,j);
    diag_term[2] += Doo(l,j);

  }
  if (delta_jl && delta_km) {
    //term0
    diag_term[1] -= Doo(l,i);
    diag_term[2] += Doo(l,i);
  }

  
  part[0].add_diag(diag_term[0]);
  part[1].add_diag(diag_term[1]);
  part[2].add_diag(diag_term[2]);


  
  return true;
}

// Chem. Phis. 329, p. 19, Table A.6
bool Zeroth_order_cap_triplet::block_ijk_lmn(const Config &row, const Config &col, Blas_matrix &block)
{
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  unsigned k = row.occ[2];
  
  unsigned l = col.occ[0];
  unsigned m = col.occ[1];
  unsigned n = col.occ[2];
  
  bool delta_il = i == l;
  bool delta_jm = j == m;
  bool delta_kn = k == n;

  bool delta_jl = j == l;
  bool delta_km = k == m;
  bool delta_jn = j == n;

  bool delta_im = i == m;
  bool delta_sym = sym(row.vir) == sym(col.vir);
  

  if (!delta_sym) return false;

  if (
      !(delta_kn && delta_jm && delta_il) &&
      !(delta_kn && delta_jm) &&
      !(delta_kn && delta_jl) &&
      !(delta_kn && delta_il) &&
      !(delta_kn && delta_im) &&
      !(delta_jm && delta_il) &&
      !(delta_jn && delta_im) &&
      !(delta_km && delta_jl) &&
      !(delta_km && delta_il) &&
      !(delta_jn && delta_il) 
      )
    return false;

  

  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(3 * row_block_size, 3 * col_block_size);
  block = 0.;
  
  Submatrix part[][3] = {
    {block(0, row_block_size, 
	   0, col_block_size),
     block(0, row_block_size, 
	   col_block_size, col_block_size),
     block(0, row_block_size, 
	   2 * col_block_size, col_block_size)},

    {block(row_block_size, row_block_size, 
	   0, col_block_size),
     block(row_block_size, row_block_size, 
	   col_block_size, col_block_size),
     block(row_block_size, row_block_size, 
	   2 * col_block_size, col_block_size)},

    {block(2 * row_block_size, row_block_size, 
	   0, col_block_size),
     block(2 * row_block_size, row_block_size, 
	   col_block_size, col_block_size),
     block(2 * row_block_size, row_block_size, 
	   2 * col_block_size, col_block_size)},
  };


  double diag_term[][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

  if (delta_kn && delta_jm && delta_il) {

    part[0][0].daxpy(1., Dvv[sym_vir]);
    part[1][1].daxpy(1., Dvv[sym_vir]);
    part[2][2].daxpy(1., Dvv[sym_vir]);

    diag_term[0][0] += gs_expt_value_;
    diag_term[1][1] += gs_expt_value_;
    diag_term[2][2] += gs_expt_value_;
  }

  if (delta_kn && delta_jm) {
    double d_li = Doo(l,i);
    diag_term[0][0] -= d_li;
    //       [1][0]
    //       [2][0]

    //       [0][1]
    diag_term[1][1] -= d_li;
    //       [2][1]

    //       [0][2]
    //       [1][2]
    diag_term[2][2] -= d_li;

  }
   
  if (delta_kn && delta_jl) {
    double d_mi = Doo(m,i);
    diag_term[0][0] += d_mi;
    //       [1][0]
    //       [2][0]
    
    //       [0][1]
    //       [1][1]
    diag_term[2][1] += d_mi;
    
    //       [0][2]
    diag_term[1][2] += d_mi;
    //       [2][2]

  }
  if (delta_kn && delta_il) {
    double d_mj = Doo(m,j);
    diag_term[0][0] -= d_mj;
    //       [1][0]
    //       [2][0]
    
    //       [0][1]
    diag_term[1][1] -= d_mj;
    //       [2][1]

    //       [0][2]
    //       [1][2]
    diag_term[2][2] -= d_mj;

  }
  if (delta_kn && delta_im) {
    double d_lj = Doo(l,j);
    diag_term[0][0] += d_lj;
    //       [1][0]
    //       [2][0]

    //       [0][1]
    //       [1][1]
    diag_term[2][1] += d_lj;

    //       [0][2]
    diag_term[1][2] += d_lj;
    //       [2][2]

  }
  if (delta_jm && delta_il) {
    double d_nk = Doo(n,k);
    diag_term[0][0] -= d_nk;
    //       [1][0]
    //       [2][0]

    //       [0][1]
    diag_term[1][1] -= d_nk;
    //       [2][1]

    //       [0][2]
    //       [1][2]
    diag_term[2][2] -= d_nk;

  }
  if (delta_jn && delta_im) {
    double d_lk = Doo(l,k);
    //       [0][0]
    diag_term[1][0] -= d_lk;
    //       [2][0]

    //       [0][1]
    //       [1][1]
    diag_term[2][1] -= d_lk;

    diag_term[0][2] -= d_lk;
    //       [1][2]
    //       [2][2]
  }
  if (delta_km && delta_jl) {
    double d_ni = Doo(n,i);
    //       [0][0]
    //       [1][0]
    diag_term[2][0] -= d_ni;

    diag_term[0][1] -= d_ni;
    //       [1][1]
    //       [2][1]

    //       [0][2]
    diag_term[1][2] -= d_ni;
    //       [2][2]
  }
  if (delta_km && delta_il) {
    double d_nj = Doo(n,j);
    //       [0][0]
    diag_term[1][0] += d_nj;
    //       [2][0]

    diag_term[0][1] += d_nj;
    //       [1][1]
    //       [2][1]

    //       [0][2]
    //       [1][2]
    diag_term[2][2] += d_nj;

  }
  if (delta_jn && delta_il) {
    double d_mk = Doo(m,k);
    //       [0][0]
    diag_term[1][0] += d_mk;
    //       [2][0]

    diag_term[0][1] += d_mk;
    //       [1][1]
    //       [2][1]

    //       [0][2]
    //       [1][2]
    diag_term[2][2] += d_mk;
    
  }
  
  

  
  
  part[0][0].add_diag(diag_term[0][0]);  part[0][1].add_diag(diag_term[0][1]); part[0][2].add_diag(diag_term[0][2]); 
  part[1][0].add_diag(diag_term[1][0]);  part[1][1].add_diag(diag_term[1][1]); part[1][2].add_diag(diag_term[1][2]);
  part[2][0].add_diag(diag_term[2][0]);  part[2][1].add_diag(diag_term[2][1]); part[2][2].add_diag(diag_term[2][2]);
  

  
  return true;  
}

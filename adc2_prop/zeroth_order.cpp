#include "zeroth_order.hpp"
#include "adc2_dip/config.hpp"
#include "matrices.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

static const double SQRT_2 =  1.41421356237310;
static const double SQRT_1_2 = 0.707106781186548;
static const double SQRT_3_2 = 1.22474487139159;
static const double SQRT_3_4 = 0.866025403784439;


// JUST A TEST!

Zeroth_order_cap::Zeroth_order_cap(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat) 
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
	  //= 100 * vir_row + vir_col;
	  Dvv[sym](sym_block_location_[vir_row], sym_block_location_[vir_col]) 
	    =  (*mat)(vir_row,vir_col);
    }
  }



//   for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) {
//     cout << "Dvv[" << sym << "]:\n";
//     Dvv[sym].print();
//   }
    
//   for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) {
//     cout << "Dvo[" << sym << "]:\n";
//     Dvo[sym].print();
//   }
  
//   cout << "Doo:\n";
//   Doo.print();
  

  

}


bool Zeroth_order_cap::block_ii_jj(const Config &row, const Config &col, double &element)
{

  
  unsigned i = row.occ[0];
  unsigned j = col.occ[0];

  bool delta_ij = i == j;
  
  element = 0.;

  if(delta_ij) 
    element += gs_expt_value_- 2. * Doo(i,i);

  return true;
  
}


bool Zeroth_order_cap::block_ij_kk(const Config &row, const Config &col, double &element)
{

  unsigned i = row.occ[0];
  unsigned j = row.occ[1];  
  unsigned k = col.occ[0];

  bool delta_jk = j == k;
  bool delta_ik = i == k;

   
  element = 0.;
  
  if (delta_jk) {
    element -= SQRT_2 * Doo(k,i);
  }
  if (delta_ik) {
    element -= SQRT_2 * Doo(k,j);
  }

  return true;
  
}

bool Zeroth_order_cap::block_ij_kl(const Config &row, const Config &col, double &element)
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
    element -= Doo(l,i);
  if (delta_il)
    element -= Doo(k,j);

  return true;
  
}


bool Zeroth_order_cap::block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned int l = row.occ[0];
  unsigned int k = row.occ[1];
  unsigned int i = col.occ[0];
  
  bool delta_ik = i == k;
  
  if (!delta_ik) return false;
  
  unsigned int block_size = size_vir_group(row.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(block_size);

  block = 0.;
  Submatrix submat = block(0,block_size, 0, 1); 
  submat.daxpy(-SQRT_2, Dvo[sym_vir](0, block_size, l, 1));

  
  return true;
}

bool Zeroth_order_cap::block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block)
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
      submat.daxpy(1., Dvo[sym_vir](0, block_size, k, 1));  

  return true;
}


bool Zeroth_order_cap::block_klm_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  return false;
}


bool Zeroth_order_cap::block_klm_ij(const Config &row, const Config &col, Blas_matrix& block)
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
  block.allocate(2 * block_size);
  block = 0.;

  Submatrix part[] = {block(0, block_size, 0, 1),
		      block(block_size, block_size, 0, 1)};
  if (delta_jm && delta_il) {
    
    part[0].daxpy(SQRT_1_2, Dvo[sym_vir](0, block_size, k, 1));
    part[1].daxpy(SQRT_3_2, Dvo[sym_vir](0, block_size, k, 1));

  }

  if (delta_jm && delta_ik) 
    part[0].daxpy(-2. * SQRT_1_2, Dvo[sym_vir](0, block_size, l, 1));

  
  
  if (delta_jl && delta_ik) {
    part[0].daxpy(SQRT_1_2, Dvo[sym_vir](0, block_size, m, 1));
    part[1].daxpy(-SQRT_3_2, Dvo[sym_vir](0, block_size, m, 1));


  }


  return true;
}


bool Zeroth_order_cap::block_jii_lkk(const Config &row, const Config &col, Blas_matrix& block)
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
      !(delta_jk && delta_il) &&
      !(delta_ik)             
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

  
  if (delta_jk && delta_il) {
    diag_term +=  Doo(k,i);
  }
  
  if (delta_ik) 
    diag_term -= Doo(l,j);

  block.add_diag(diag_term);

  
  return true;
}


bool Zeroth_order_cap::block_ijk_mll(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  unsigned k = row.occ[2];
  
  unsigned m = col.occ[0];
  unsigned l = col.occ[1];
  

  bool delta_jm = j == m;
  bool delta_il = i == l;
  bool delta_km = k == m;
  bool delta_im = i == m;
  bool delta_kl = k == l;
  bool delta_jl = j == l;

  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  if (!delta_sym) return false;

  if (
      !(delta_jm && delta_kl) &&
      !(delta_im && delta_jl) &&
      !(delta_im && delta_kl) &&
      !(delta_jl && delta_km) &&
      !(delta_il && delta_jm) && 
      !(delta_il && delta_km) 
      )
    return false;

  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  unsigned sym_vir = sym(row.vir);
  block.allocate(2 * row_block_size, col_block_size);
  block = 0.;
  
  Submatrix part[] = 
    {block(0, row_block_size, 0, col_block_size),
     block(row_block_size, row_block_size, 0, col_block_size)};

  double diag_term[] = {0., 0.};
  
  if (delta_jm && delta_kl) {
    diag_term[0] -= 2. * SQRT_1_2 * Doo(l,i);
    //diag_term[1] 
  }
  if (delta_im && delta_jl) {
    diag_term[0] += SQRT_1_2 * Doo(l,k);
    diag_term[1] += SQRT_3_2 * Doo(l,k);
  }
  if (delta_im && delta_kl) {
    diag_term[0] += SQRT_1_2 * Doo(l,j);
    diag_term[1] += SQRT_3_2 * Doo(l,j);
  }
  //Aditional:
  if (delta_jl && delta_km) {
    diag_term[0] += SQRT_1_2 * Doo(l,i);
    diag_term[1] -= SQRT_3_2 * Doo(l,i);    
  }
  if (delta_il && delta_jm) {
    diag_term[0] -= SQRT_1_2 * 2.* Doo(l,k);
    //[1]
  }
  if (delta_il && delta_km) {
    diag_term[0] += SQRT_1_2 * Doo(l,j);
    diag_term[1] -= SQRT_3_2 * Doo(l,j);
  }


  part[0].add_diag(diag_term[0]);
  part[1].add_diag(diag_term[1]);
  
  return true;
}

// Chem. Phis. 329, p. 19, Table A.6
bool Zeroth_order_cap::block_ijk_lmn(const Config &row, const Config &col, Blas_matrix &block)
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
  block.allocate(2 * row_block_size, 2 * col_block_size);
  block = 0.;
  
  Submatrix part[][2] = {
    {block(0, row_block_size, 0, col_block_size),
     block(0, row_block_size, col_block_size, col_block_size)},
    
    {block(row_block_size, row_block_size, 0, col_block_size),
     block(row_block_size, row_block_size, col_block_size, col_block_size)}};


  double diag_term[][2] = {{0., 0.}, {0., 0.}};

  if (delta_kn && delta_jm && delta_il) {

    part[0][0].daxpy(1., Dvv[sym_vir]);
    part[1][1].daxpy(1., Dvv[sym_vir]);

    diag_term[0][0] += gs_expt_value_;
    diag_term[1][1] += gs_expt_value_;
  }

  if (delta_kn && delta_jm) {
    diag_term[0][0] -= Doo(l,i);
    diag_term[1][1] -= Doo(l,i);

  }
    
  if (delta_kn && delta_jl) {
    diag_term[0][0] += 0.5 * Doo(m,i);
    diag_term[0][1] += SQRT_3_4 * Doo(m,i);
    diag_term[1][0] += SQRT_3_4 * Doo(m,i);
    diag_term[1][1] -= 0.5 * Doo(m,i);

  }
  if (delta_kn && delta_il) {
    diag_term[0][0] -= Doo(m,j);
    diag_term[1][1] -= Doo(m,j);

  }
  if (delta_kn && delta_im) {
    diag_term[0][0] += 0.5 * Doo(l,j);
    diag_term[0][1] += SQRT_3_4 * Doo(l,j);
    diag_term[1][0] += SQRT_3_4 * Doo(l,j);
    diag_term[1][1] -= 0.5 * Doo(l,j);

  }
  if (delta_jm && delta_il) {
    diag_term[0][0] -= Doo(n,k);
    diag_term[1][1] -= Doo(n,k);

  }
  if (delta_jn && delta_im) {
    diag_term[0][0] += 0.5 * Doo(l,k);
    diag_term[0][1] -=  SQRT_3_4 * Doo(l,k);
    diag_term[1][0] +=  SQRT_3_4 * Doo(l,k);
    diag_term[1][1] += 0.5 * Doo(l,k);

  }
  if (delta_km && delta_jl) {
    diag_term[0][0] += 0.5 * Doo(n,i);
    diag_term[0][1] += SQRT_3_4 * Doo(n,i);
    diag_term[1][0] -= SQRT_3_4 * Doo(n,i);
    diag_term[1][1] += 0.5 * Doo(n,i);

  }
  if (delta_km && delta_il) {
    diag_term[0][0] += 0.5 * Doo(n,j);
    diag_term[0][1] -= SQRT_3_4 * Doo(n,j);
    diag_term[1][0] -= SQRT_3_4 * Doo(n,j);
    diag_term[1][1] -= 0.5 * Doo(n,j);
 
  }
  if (delta_jn && delta_il) {
    diag_term[0][0] += 0.5 * Doo(m,k);
    diag_term[0][1] -= SQRT_3_4 * Doo(m,k);
    diag_term[1][0] -= SQRT_3_4 * Doo(m,k);
    diag_term[1][1] -= 0.5 * Doo(m,k);
  }
  
  
  
  part[0][0].add_diag(diag_term[0][0]);  part[0][1].add_diag(diag_term[0][1]);
  part[1][0].add_diag(diag_term[1][0]);  part[1][1].add_diag(diag_term[1][1]);
  

  
  return true;  
}

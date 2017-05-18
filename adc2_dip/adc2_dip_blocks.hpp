#ifndef __ADC2_DIP_BLOCKS_HPP__
#define __ADC2_DIP_BLOCKS_HPP__


// The file contains the declaration of the ADC2_DIP_blocks class

// The abstract ADC2_DIP_blocks class provides an interface for building the ADC2 matrix blocks.
// It also contains the static SCF data (and protected methods operating with it) used in the construction of the blocks
// shared by both triplet and singlet. For a definition of Vij,k, Aij and Bij integral matrices see
// Tarantelli's article, page 13.
// 
// The implementation of the concrete singlet and triplet formulas is done by its derived classes.


#include "scf_data/scf_data_reader.hpp"
#include "blas_matrix.hpp"
#include "integral_table.hpp"



struct Config;


class ADC2_DIP_blocks {    

  SCF_data_reader& phis_;                // Provides the SCF data

  static bool is_static_set;

  static std::vector<unsigned int> vir_group_sizes_;  // Stores the number of virtual orbitals in a given symmetry.
                                                      // Basically, the size of the building blocks is a multiple of one of these numbers
  static std::vector<Blas_matrix> diag_energies_; // Stores the virtual-orbital energies, sorted by symmetry,
                                                  // that are to be added to the diagonals of some of the building blocks
  static std::vector<double> orb_energy_;         // Stores the orbital energies

protected:   // The protected section contains methods for accessing the static SCF data
  static Integral_table* int_block_;     // Provides the sorted two-electron integrals

  inline static double V1122(unsigned a, unsigned b, unsigned c, unsigned d) 
  { return int_block_->integral(a, b, c, d); } // get a two-electron integral

  inline static double V1122_MINUS(unsigned a, unsigned b, unsigned c, unsigned d)
  { return V1122(a, b, c, d) - V1122(a, d, b, c); }
  
  inline static double V1122_PLUS(unsigned a, unsigned b, unsigned c, unsigned d) 
  { return V1122(a, b, c, d) + V1122(a, d, b, c); }

  inline unsigned sym(unsigned orb)  
  { return phis_.irrep(orb); } // get the symmetry of an orbital

  inline unsigned sym_product(unsigned sym1, unsigned sym2)
  { return phis_.irrep_product(sym1,sym2);        } // Get the product of two symmetries

  inline unsigned sym_product(unsigned sym1, unsigned sym2, unsigned sym3)
  { return sym_product(sym_product(sym1,sym2),sym3); } // ... of three symmetries

  inline unsigned sym_product(unsigned sym1, unsigned sym2, unsigned sym3, unsigned sym4)
  { return sym_product(sym_product(sym_product(sym1,sym2),sym3), sym4); }  // ... of four symmetries

  inline unsigned number_occupied() 
  { return phis_.number_occupied(); } // Get the number of occupied orbitals

  inline unsigned number_orbitals()
  { return phis_.number_orbitals(); } // Get the total number of orbitals

  inline static double energy(unsigned orb)
  { return orb_energy_[orb]; } // Get the orbital energy

  inline static Daxpy_argument* A(unsigned i, unsigned j, unsigned sym)
  { return int_block_->get_risj_block(i,j,sym); } // Get a (ri|sj) block

  inline static Daxpy_argument* B(unsigned i, unsigned j, unsigned sym)
  { return int_block_->get_rsij_block(i,j,sym); } // Get a (rs|ij) block
  
  inline static Daxpy_argument* V(unsigned i, unsigned j, unsigned k)
  { return int_block_->get_rkij_block(i,j,k); }   // Get a (rk|ij) block

  inline static Blas_matrix& diag_energies(unsigned sym)
  { return diag_energies_[sym];  } // Get the diagonal term of virtual energies

public:

  ADC2_DIP_blocks(SCF_data_reader& phis, Integral_table& tab);

  
  // The functions responsible for generating the building block 
  // of the ADC2 matrix are overridden by singlet and triplet methods
  virtual bool block_ii_jj(const Config& row, const Config& col, 
			   double &element) = 0;
  virtual bool block_ij_kk(const Config& row, const Config& col, 
			   double &element) = 0;
  virtual bool block_ij_kl(const Config& row, const Config& col, 
			   double &element) = 0;
  virtual bool block_lkk_ii(const Config& row, const Config& col, 
			    Blas_matrix& block) = 0;
  virtual bool block_lkk_ij(const Config& row, const Config& col, 
			    Blas_matrix& block) = 0;
  virtual bool block_klm_ii(const Config& row, const Config& col, 
			    Blas_matrix& block) = 0;
  virtual bool block_klm_ij(const Config& row, const Config& col, 
			    Blas_matrix& block) = 0;
  virtual bool block_jii_lkk(const Config& row, const Config& col, 
			     Blas_matrix& block) = 0;
  virtual bool block_ijk_mll(const Config& row, const Config& col, 
			     Blas_matrix& block) = 0;
  virtual bool block_ijk_lmn(const Config& row, const Config& col, 
			     Blas_matrix& block) = 0;  

  inline unsigned size_vir_group(unsigned vir) 
  {return vir_group_sizes_[sym(vir)];}  // Get the number of virtual orbitals for a given symmetry
};


#endif //#ifndef __ADC2_DIP_BLOCKS_HPP__

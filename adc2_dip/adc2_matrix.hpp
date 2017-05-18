#ifndef __ADC2_MATRIX_HPP__
#define __ADC2_MATRIX_HPP__

// The ADC2_matrix class represents the ADC2 double ionization matrix
// 
// Main use: to perform the matrix vector multiply step
// used in the Lanczos algorithm.
// The structure of the matrix is based on 
// F. Tarantelli, Chemical Physics 329 (2005) 11-21, section 3.2


#include "../adc_matrix.hpp"
#include "config.hpp"
#include <vector>


class ADC_analyzer;
class SCF_data_reader;
class ADC2_DIP_blocks;
class Blas_matrix;
class Integral_table;

class Adc2_matrix: public ADC_matrix {
  
  SCF_data_reader& phis_;                    // provides the scf data
  ADC2_DIP_blocks* blocks_;                  // provides the matrix building blocks
  unsigned int count_mult;                   // used to emulate subspace iteration, see adc2_matrix.cpp
  
  unsigned int symmetry_;                    // holds the symmetry of the matrix
  unsigned int multiplet_functions_;         // holds the number of |3h1p> types for singlet/triplet
                                             // see Table A.1. page 17, important for the sizes of the building blocks

  std::vector<Config> configs_;              // stores the configurations relevant to the current symmetry and spin
                                             // basically, this vector is a map between the physical indexes of the matrix
                                             // and the corresponding configurations
  
  unsigned int begin_ij_;                    // marks the beginning of the |ij> configurations 
  unsigned int begin_jii_;                   // marks the beginning of the |jiir> configurations
  unsigned int begin_ijk_;                   // marks the beginning of the |ijkr>

  std::vector<unsigned int> jii_, ijk_;      // stores the starting indexes of the groups of |3h1p> configurations with
                                             // the same hole orbitals
  
  void add_ii_configs();                     // These functions generate the electronic configurations 
  void add_ij_configs();                     // corresponding to the rows and columns of the ADC2 matrix
  void add_jiir_configs();
  void add_ijkr_configs();

  void build_submatrix_ii_jj(Blas_matrix &adc2_matrix);      
  void build_submatrix_ij_kk(Blas_matrix &adc2_matrix);       // These functions construct the full ADC2 matrix
  void build_submatrix_ij_kl(Blas_matrix &adc2_matrix);       // from the individual blocks
  void build_submatrix_lkk_ii(Blas_matrix &adc2_matrix);
  void build_submatrix_lkk_ij(Blas_matrix &adc2_matrix);
  void build_submatrix_klm_ii(Blas_matrix &adc2_matrix);
  void build_submatrix_klm_ij(Blas_matrix &adc2_matrix);
  void build_submatrix_jii_lkk(Blas_matrix &adc2_matrix);
  void build_submatrix_ijk_mll(Blas_matrix &adc2_matrix);
  void build_submatrix_ijk_lmn(Blas_matrix &adc2_matrix);
  
  void multiply_submatrix_ii_jj(Blas_matrix** block_out, Blas_matrix** block_in, int count);  // These functions perform 
  void multiply_submatrix_ij_kk(Blas_matrix** block_out, Blas_matrix** block_in, int count);  // the vector matrix multiplication with 
  void multiply_submatrix_ij_kl(Blas_matrix** block_out, Blas_matrix** block_in, int count);  // each of the subblocks of a given type without
  void multiply_submatrix_lkk_ii(Blas_matrix** block_out, Blas_matrix** block_in, int count); // building the whole matrix
  void multiply_submatrix_lkk_ij(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  void multiply_submatrix_klm_ii(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  void multiply_submatrix_klm_ij(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  void multiply_submatrix_jii_lkk(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  void multiply_submatrix_ijk_mll(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  void multiply_submatrix_ijk_lmn(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  
  void print_configs();
 
public:

  Adc2_matrix(SCF_data_reader& phis, Integral_table& tab, unsigned int sym, unsigned int spin);
  virtual ~Adc2_matrix();
  
  virtual unsigned int size() const {return configs_.size();}                          // returns the size of the matrix
  virtual int operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count);  // the matrix vector multiply routine
  virtual void build_matrix(Blas_matrix& mat);                                         // Generates the whole ADC2 matrix
  
  virtual int accept_analyzer(ADC_analyzer& analyst);                                 // Accepts a class that will perform the analysis 
                                                                                       // of the matrix's properties
  
  virtual unsigned int main_block_size() const {return begin_jii_;}                          // returns the size of the main part
  virtual std::string get_conf(unsigned int i) const;                                        // returns a string containing the configuration
                                                                                       // mapped to the given row/column of the matrix
  virtual unsigned int symmetry() const {return symmetry_;}                                  // returns the symmetry of the matrix
  virtual unsigned int spin() const {return ( multiplet_functions_ == 2 ) ? 0 : 2;}          // returns the spin of the matrix
  
  void reset(){count_mult = 0;}
};

#endif //#ifndef __ADC2_MATRIX_HPP__

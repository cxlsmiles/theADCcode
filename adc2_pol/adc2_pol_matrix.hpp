#ifndef __ADC2_POL_MATRIX_HPP__
#define __ADC2_POL_MATRIX_HPP__



#include "adc_matrix.hpp"
#include "blas_matrix.hpp"
#include <set>


class SCF_data_reader;
//class Integral_table;

typedef std::set<unsigned> ListOfOrbs ;


class ADC2_pol_matrix: public ADC_matrix {
  SCF_data_reader* phis_;       // Provides the SCF data

  int sym_;

  int main_size, sat_size;
  int* kpq;

  int count_mult;


  void multiply_main_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count);
  void multiply_main_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count) ;
  void multiply_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count);
  
public:
  
  ADC2_pol_matrix(SCF_data_reader& phis, unsigned int sym, ListOfOrbs holes);
  virtual ~ADC2_pol_matrix();
  
  virtual void build_matrix(Blas_matrix& mat);
  virtual int operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count);
  virtual unsigned int main_block_size() const {return main_size;}
  virtual std::string get_conf(unsigned int i) const;

  virtual int accept_analyzer( ADC_analyzer& adc_an );
  virtual unsigned int size() const {return main_size+sat_size;}
  virtual unsigned int symmetry() const {return sym_;}
  virtual unsigned int spin() const {return 1;}
  virtual void reset() {count_mult = 0;}
};



#endif //#ifndef __ADC2_POL_MATRIX_HPP__

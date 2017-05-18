#ifndef __ND_ADC3_CAP_MATRIX_HPP__
#define __ND_ADC3_CAP_MATRIX_HPP__


// This file declares the nd_adc3_matrix class which
// represents a (N-1) ND_ADC3 propagator.
// The structure of the matrix follows Joerg Breidbach's nd_adc3 code
// and J. Chem. Phys. 109 (1998) 4734


#include "adc_matrix.hpp"
#include "blas_matrix.hpp"
#include "matrices.hpp"
#include <vector>

struct symtab_t; // These structures are required by Joerg's code
struct scf_t;    // and contain some sorted scf information
struct inp_t;

class SCF_data_reader;
// class Integral_table;


// These functions supply Joerg's code with the two electron integrals
// and the multiplication tables contained in the nd_adc3_matrix class.
// Below they are declared as friends of the class 
// extern "C" {
//   double Vpqrs_cap(int*, int*, int*, int*);
//   int Multab_cap(int,int);

// }

class ND_ADC3_CAP_matrix: public ADC_matrix {
  static SCF_data_reader* phis_;       // Provides the SCF data
  //  static Integral_table *int_tab;      // Provides the two electron integrals.


  int dim_1h;                          // Stores the dimensions of the main space
  int dim_2h1p;                        // Stores the dimensions of the satellite states

  unsigned sym_;                       // Stores the symmetry of the propagator
  int count_mult;                      // Counts the number of Lanczos iterations see operator()

  std::vector<std::string> configs_;   // Stores electron configuration labels to each row/column


  symtab_t* symtab;   // Structures needed for calling Joerg's functions
  scf_t* scf;
  inp_t* inp;
  
  Triangular_matrix<double>* cap_mat;
  static Triangular_matrix<double>* dens_mat;
  double d_null; // Second order expectation value 
  double d_null_null; // In Imke's notation, the zeroth order expectation value of the CAP operator

  void scf_data(scf_t* scf); // Initializes the scf structure suing the SCF_data reader
  void add_configs();        // Initializes the configuration labels 

  // Methods building parts of the ND ADC3 matrix
  void build_main_block(Blas_matrix& mat);  
  void build_main_sat_block(Blas_matrix& mat);
  void build_sat_block(Blas_matrix& mat);

  //Methods multiplying parts of the ND ADC3 matrix with a set of vectors
//   void multiply_main_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count);
//   void multiply_main_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count);
//   void multiply_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count);
//   void multiply_sat_row(Blas_matrix** block_out,  Blas_matrix** block_in, int count,
// 			Submatrix row, unsigned block_rows);
public:
  
  ND_ADC3_CAP_matrix(SCF_data_reader& phis, Triangular_matrix<double>* mat, unsigned int sym, unsigned int mode);
  virtual ~ND_ADC3_CAP_matrix();

  virtual void build_matrix(Blas_matrix& mat);
  virtual int operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count) {return 0;}
  virtual unsigned int main_block_size() const {return dim_1h;}
  virtual std::string get_conf(unsigned int i) const {return configs_[i];}
  virtual int accept_analyzer( ADC_analyzer& adc_an );
  virtual unsigned int size() const {return dim_1h + dim_2h1p;}
  virtual unsigned int symmetry() const {return sym_;}
  virtual unsigned int spin() const {return 1;}
  virtual void reset(){count_mult = 0;}

//   friend double Vpqrs_cap(int* a, int* b, int* c, int* d);
//   friend int    Multab_cap(int a, int b);
};



#endif //#ifndef __ND_ADC3_CAP_MATRIX_HPP__

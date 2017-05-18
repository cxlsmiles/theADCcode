#ifndef __ADC_MATRIX_HPP__
#define __ADC_MATRIX_HPP__

#include <string>
class ADC_analyzer;
class Blas_matrix;


// This class declares the interface for
// ADC propagator's implementations.
// These must be able to build a whole matrix, or perform
// at a matrix-vector multiply that is used in the Lanczos diagonalization procedure


class ADC_matrix {

public:


  virtual ~ADC_matrix() {}
  
  // This method accepts an external visitor
  // that is able to use the properties of the ADC matrix for a specific analysis
  virtual int accept_analyzer( ADC_analyzer& adc_an ) = 0;

  // A matrix-vector multiply: Here a block of input/output vectors are used
  virtual int operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count) = 0;
  // Builds the whole ADC matrix
  virtual void build_matrix(Blas_matrix& mat) = 0;
  
  virtual unsigned int size() const = 0;      // Returns the size of the matrix
  virtual unsigned int main_block_size() const  = 0; // Returns the size of its main space
  virtual unsigned int symmetry() const = 0;        // Returns the symmetry for which the matrix was build 
  virtual unsigned int spin() const = 0;            // Returns the spin

  virtual std::string get_conf(unsigned int i) const = 0; // Returns the corresponding electronic configuration (a string)
                                                    // for a given row/column of the ADC matrix
  virtual void reset() = 0;
};

#endif //#ifndef __ADC_MATRIX_HPP__

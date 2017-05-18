#ifndef __ADC_ANALYZER_HPP__
#define __ADC_ANALYZER_HPP__

// This file contains the declaration of a class that
// analyzes the properties of the ADC propagators. 
// It produces the spectrum of a given propagator. 
// Its derived classes can define a different type of analysis

#include <vector>


class ADC_matrix;
class Blas_matrix;
struct Diag_info;
struct Eigen_info;



class ADC_analyzer {
  
  struct Diag_info& diag_;
  struct Eigen_info& eig_;

  void find_spurious(int, Blas_matrix&, std::vector<unsigned>&);

protected:

  // Prints the results
  void print_eigenvals(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms);
  // Uses lapack to diagonalize the matrix
  void lapack_diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals);
  // Uses Lanczos to diagonalize the matrix
  void lanczos_diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms);
  void diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms);
 
public:
  
  ADC_analyzer(struct Diag_info& diag, struct Eigen_info& eig) : diag_(diag), eig_(eig) {}
  virtual int analyze(ADC_matrix& adc_mat);

};


#endif // #ifndef __ADC_ANALYZER_HPP__



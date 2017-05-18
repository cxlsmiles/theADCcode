#ifndef __RSCAP_ANALYZER_HPP__
#define __RSCAP_ANALYZER_HPP__

#include "adc_analyzer.hpp"

class SCF_data_reader;
class ADC_matrix;

struct Diag_info;
struct Eigen_info;
struct Cap_info;


class RSCAP_analyzer: public ADC_analyzer {
  
  SCF_data_reader& phis_;
  
  struct Cap_info& cap_;
  
  void print_nicolas_input_ip( ADC_matrix& adc_mat);
  void print_nicolas_input_dip( ADC_matrix& adc_mat);

public:

  RSCAP_analyzer(SCF_data_reader& phis, struct Diag_info& diag, 
		 struct Eigen_info& eig, struct Cap_info& cap);
  // Get the matrix that is to be analyzed
  virtual int analyze(ADC_matrix& adc_mat);

};




#endif // #ifndef __RSCAP_ANALYZER_HPP__


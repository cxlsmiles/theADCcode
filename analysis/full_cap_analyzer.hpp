#ifndef __FULL_CAP_ANALYZER_HPP__
#define __FULL_CAP_ANALYZER_HPP__

#include "adc_analyzer.hpp"


class SCF_data_reader;
class Integral_table;

struct Diag_info;
struct Eigen_info;
struct Cap_info;

class Full_CAP_analyzer: public ADC_analyzer {
  
  struct Cap_info& cap_;
  SCF_data_reader& phis_;
  Integral_table& tab_;
  
public:

  Full_CAP_analyzer(SCF_data_reader& phis, Integral_table& tab, 
		    struct Diag_info& diag, struct Eigen_info& eig, struct Cap_info& cap);
  virtual int analyze(ADC_matrix& adc_mat);
  
};




#endif // #ifndef __FULL_CAP_ANALYZER_HPP__


#ifndef __ISR_DIPOLE_ANALYZER_HPP__
#define __ISR_DIPOLE_ANALYZER_HPP__


#include "adc_analyzer.hpp"

class SCF_data_reader;
class Integral_table;

struct Diag_info;
struct Eigen_info;
struct Dip_info;


class ISR_dipole_analyzer: public ADC_analyzer {
  
  struct Dip_info& dip_;
  SCF_data_reader& phis_;
  Integral_table& tab_;
  
public:

  ISR_dipole_analyzer(SCF_data_reader& phis, Integral_table& tab, struct Diag_info& diag, 
		      struct Eigen_info& eig, struct Dip_info& dip) 
    : phis_(phis), tab_(tab), ADC_analyzer(diag, eig), dip_(dip) {}
  
  // Get the matrix that is to be analyzed
  virtual int analyze(ADC_matrix& adc_mat);
  
};




#endif // #ifndef __ISR_DIPOLE_ANALYZER_HPP__


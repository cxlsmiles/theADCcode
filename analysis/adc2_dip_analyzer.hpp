#ifndef __ADC2_DIP_ANALYZER_HPP__
#define __ADC2_DIP_ANALYZER_HPP__

#include "adc_analyzer.hpp"

#include <list>
#include <utility>
#include <set>
#include <string>


class ADC_matrix;
class Blas_matrix;
class SCF_data_reader;

struct Diag_info;
struct Eigen_info;


typedef std::list<std::pair<std::string, std::set<unsigned> > > OrbSet;

class ADC2_DIP_analyzer: public ADC_analyzer {
  
  SCF_data_reader& _phis;
  OrbSet _groups;

  void twohole_mulliken(ADC_matrix& adc2_mat, Blas_matrix& eig_vecs, Blas_matrix& eig_vals);

public:
  
  ADC2_DIP_analyzer(SCF_data_reader& phis,  OrbSet groups, 
		    struct Diag_info& diag, struct Eigen_info& eig);
  virtual int analyze(ADC_matrix& adc_mat );
  
};



#endif // #ifndef __ADC2_DIP_ANALYZER_HPP__


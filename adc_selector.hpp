#ifndef __ADC_SELECTOR_HPP__
#define __ADC_SELECTOR_HPP__


// Choose the concrete functionality according ot
// the user's input


#include <set>

typedef std::set<unsigned> ListOfSyms ;
typedef std::set<unsigned> ListOfSpins ;

class ADC_matrix;
class ADC_analyzer;
struct Input_data;
class SCF_data_reader;

class Adc_selector {

  Input_data* input;
  SCF_data_reader* reader_;

  void check_input();
  void print_input();

public:

  Adc_selector();
  ~Adc_selector();
  ADC_matrix* new_adc_matrix(unsigned sym, unsigned spin);
  ADC_analyzer* new_adc_analyzer();

  
  ListOfSyms list_of_syms();
  ListOfSpins list_of_spins();
};







#endif //#ifndef __ADC_SELECTOR_HPP__

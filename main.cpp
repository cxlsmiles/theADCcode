#include "analysis/adc_analyzer.hpp"
#include "adc_matrix.hpp"
#include "adc_selector.hpp"

#include <iostream>

int main()
{
  std::cout << "  TheADCcode: A collection of ADC/ISR source codes.\n"
       << "  Developed for the TC group, University of Heidelberg.\n"
       << "  Contributors: Yasen Velkov (yasen.velkov@pci.uni-heidelberg.de),\n"
       << "                Anthony Dutoi, Nicolas Sisourat, Tsveta Miteva,\n"
       << "                Joerg Breidbach, Imke Mueller, Nayana Vaval,\n"
       << "                Francesco Tarantelli, Soeren Kopelke,\n"
       << "                Sajeev Yesodharan, Kirill Gokhberg, Robin Santra\n\n";
 
  try {
    
    Adc_selector* selector;
    selector = new Adc_selector; // Read the input information 
    
    ADC_analyzer* diagonalizer = selector->new_adc_analyzer();
    
    ListOfSpins spins = selector->list_of_spins();
    ListOfSyms symms = selector->list_of_syms();
    
    ListOfSpins::iterator spin; 
    ListOfSyms::iterator sym;
    for(spin = spins.begin(); spin != spins.end(); spin++)     
      for(sym = symms.begin(); sym != symms.end(); sym++) {
	
	ADC_matrix* mat = selector->new_adc_matrix(*sym,*spin);
	
	mat->accept_analyzer(*diagonalizer);
	
	delete mat;
      }
    
    
    delete diagonalizer;
    delete selector;
    
  } catch (std::string er) {
    std::cout << ' ' << er << std::endl;
    std::cout << " Task not accomplished.\n";
    return 1;
  }

  std::cout << "\n Task accomplished.\n";
  return 0;

}

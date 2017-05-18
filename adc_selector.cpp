#include "adc_selector.hpp"
#include "adc2_dip/adc2_matrix.hpp"
#include "adc2_pol/adc2_pol_matrix.hpp"
#include "ndadc3_ip/nd_adc3_matrix.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "analysis/adc2_dip_analyzer.hpp"
#include "analysis/full_cap_analyzer.hpp"
#include "analysis/isr_dipole_analyzer.hpp"
#include "analysis/subspacecap_analyzer.hpp"


#include "input_data.hpp"
#include "input_struct.hpp"
#include "integral_table.hpp"


#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>



#define OOOO 1
#define OOVO 2
#define VOVO 4
#define OOVV 8
#define VOVV 16
#define VVVV 32 



using namespace std;



extern Integral_table* integral_table = 0;

Adc_selector::Adc_selector() : input(0), reader_(0)
{
    input = new Input_data();

    if (input->backend.empty()) throw string("Please specify the front-end SCF program.\n");
    reader_ = new SCF_data_reader(input->backend.c_str());

    check_input();

    int analyzer = 0;

    if (input->cap) analyzer++;
    if (input->dip) analyzer++;
    if (!input->el_groups.empty()) analyzer++;//popana
    if (analyzer > 1) throw string("Can't make more than one type of analysis... for now.\nAvailable:cap,dip,popana\n");

    //The double ionization ADC2 propagator needs only OOOO,OOVO,VOVO,OOVV type of integrals
    // Actually the corresponding N-2 ISR also needs the same type of integrals, however
    // it uses the Sigma4+ to obtain the density corrections, which on the other hand uses all types of
    // integrals. Therefore, load all integrals unless just the ADC2 DIP propagator is invoked.
    if ((input->propag->method == ADC2DIP || input->propag->method == ADC2IPx) &&  !input->cap && !input->dip )
      integral_table = new Integral_table(*reader_, OOOO|OOVO|VOVO|OOVV);
    else
      integral_table = new Integral_table(*reader_, OOOO|OOVO|VOVO|OOVV|VOVV|VVVV);

    print_input();
    
}


Adc_selector::~Adc_selector() 
{
  delete input;
}



ListOfSyms Adc_selector::list_of_syms() 
{
  return input->propag->symms;
}


ListOfSpins Adc_selector::list_of_spins() 
{
   return input->propag->spin;
}


ADC_analyzer* Adc_selector::new_adc_analyzer()
{
  if (input->cap) {
    
    if (input->cap->type == SUBSPACE) {
      return new RSCAP_analyzer(*reader_, *(input->diag), *(input->eig), *(input->cap));
    }else if (input->cap->type == FULL){
      return new Full_CAP_analyzer(*reader_, *integral_table, *(input->diag), *(input->eig), *(input->cap));
    } else 
      throw string ("Unknown CAP type.\n");
  } else if (input->dip) { 
    return new ISR_dipole_analyzer(*reader_, *integral_table, *(input->diag), *(input->eig), *(input->dip));
  }else if (!input->el_groups.empty())
    return new ADC2_DIP_analyzer(*reader_, input->el_groups, *(input->diag), *(input->eig));
  else
    return new ADC_analyzer(*(input->diag), *(input->eig));
}


ADC_matrix* Adc_selector::new_adc_matrix(unsigned int sym, unsigned int spin)
{
  switch (input->propag->method) {
  case ADC2DIP:
    return new Adc2_matrix(*reader_, *integral_table, sym, spin);
  case ADC2PP:
    return new ADC2_pol_matrix(*reader_, sym, input->propag->holes);    
  case ADC2IPx:
    return new ND_ADC3_matrix(*reader_, sym, spin, 2);
  case NDADC3IP:
    return new ND_ADC3_matrix(*reader_, sym, spin, 0);
  case NDADC3AP:
    return new ND_ADC3_matrix(*reader_, sym, spin, 1);
  default:
    throw string("Unknown propagator selected.\n");
  }
}


// Check the consistency of the input
void Adc_selector::check_input() 
{
   
  if(!input->propag) throw string("Propagator input not specified.\n");
  if (input->propag->spin.empty()) throw string("Spin input not specified.\n");

  // Default settings
  if (!input->diag) { input->diag = new struct Diag_info; input->diag->type=DiagFull;}
  if (!input->eig) { input->eig = new struct Eigen_info; input->eig->ps = 0.01; input->eig->thresh = 0.0001;}
  set<unsigned>::iterator it;
  if (input->propag->symms.empty()) {
    for (int i = 0; i < reader_->number_irreps(); i++)
      input->propag->symms.insert(i); 
  }
}

void Adc_selector::print_input() 
{

  cout << " Input summary:\n";
  

  string method;
  switch (input->propag->method) {
  case ADC2DIP: method.assign(" Second-order ADC propagator for double ionization potentials"); break;
  case ADC2PP: method.assign(" Second-order ADC polarization propagator"); break;
  case ADC2IPx: method.assign(" Second-order extended ADC propagator for ionization potentials"); break;
  case NDADC3IP: method.assign(" Third-order Non-Dyson ADC propagator for ionization potentials"); break;
  case NDADC3AP: method.assign(" Third-order Non-Dyson ADC propagator for affinity potentials"); break;
  default: method.assign(" Unknown method");
  }
  cout << method << endl;

  cout << setw(25) << " Symmetries:";
  set<unsigned>::iterator it;
  for (it = input->propag->symms.begin(); it != input->propag->symms.end(); it++)
    cout << ' ' << *it+1;
  cout << endl;

  cout << setw(25)  << " Spins:";
  for (it = input->propag->spin.begin(); it != input->propag->spin.end(); it++)
    cout << ' ' << *it+1;
  cout << endl;
  
  cout << setw(25)  << " Diagonalization:";
  switch (input->diag->type) {
  case DiagFull: cout << " Full (LAPACK)"; break;
  case DiagLanc: cout << " Lanczos"; break;
  default: cout << " Unknown";
  }
  cout << endl;


  bool task_printed = false;
  cout  << setw(25) << "Task:";
  if (input->cap) {
    cout << ' ' << 
      ((input->cap->type == FULL) ? "Full" : "Subspace") << "-CAP trajectories\n";
    task_printed = true;
  }
  if (!input->el_groups.empty())  {
    cout <<  " Two-hole population analysis\n";
    task_printed = true;
  }
  if (input->dip) {
    cout  << " Compute Dipole moments (TEST ONLY)\n";
    task_printed = true;
  }
  if (!task_printed) cout << " Print poles\n";

}

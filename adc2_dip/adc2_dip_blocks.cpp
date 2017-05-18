#include "adc2_dip_blocks.hpp"

// These macros define the integrals that are 
// to be loaded by the integral table

// Define the static members  
Integral_table* ADC2_DIP_blocks::int_block_ = 0;

bool ADC2_DIP_blocks::is_static_set = false;

std::vector<unsigned int> ADC2_DIP_blocks::vir_group_sizes_;
std::vector<Blas_matrix> ADC2_DIP_blocks::diag_energies_;  
std::vector<double> ADC2_DIP_blocks::orb_energy_;



// the constructor initialized the static SCF data the first time it runs
ADC2_DIP_blocks::ADC2_DIP_blocks(SCF_data_reader& phis, Integral_table& tab): phis_(phis)
{

  int_block_ = &tab;
  if(!is_static_set) {
    
    //stores the position of the virtual orbitals for each symmetry group
    std::vector<unsigned int> sym_block_location_;
   
    vir_group_sizes_.assign(phis_.number_irreps(), 0);
    sym_block_location_.assign(phis_.number_orbitals(),0);
    
    // Find the size of each symmetry group of virtual orbitals
    // and map the virtual orbitals to their position in the group
    for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals(); i++) 
      sym_block_location_[i] = vir_group_sizes_[phis_.irrep(i)]++;
    
    // Set the diagonal terms of virtual energies
    for(unsigned int sym = 0; sym < phis_.number_irreps(); sym++) {
      diag_energies_.push_back(Blas_matrix(vir_group_sizes_[sym]));
      for(unsigned int vir = phis_.number_occupied(); vir < phis_.number_orbitals(); vir++)
	if (phis_.irrep(vir) == sym)
	  diag_energies_[sym](sym_block_location_[vir]) 
	    = phis_.energy(vir);
    }

    // Get the energies of all orbitals
    for(unsigned i = 0; i < phis_.number_orbitals(); i++) 
      orb_energy_.push_back(phis_.energy(i));
    
    is_static_set = true;
    
  }
  
}





#ifndef __SELF_ENERGY_HPP__
#define __SELF_ENERGY_HPP__


#include <vector>


// This class provides an interface and implementations
// for computing the static self energy (\Sigma(\inf)) and 
// the one-particle density matrix \rho. 
// Refs.: 
// J. Schirmer, A.B. Trofimov, and G.Stelter, J. Chem. Phis 109, 4734 (1998)
// A.B. Trofimov, J. Schirmer, J. Chem. Phis. 123, 144115 (2005)
// W. von Niessen et al., Computer Physics Reports 1, 57-125 (1984).


class Integral_table;
class SCF_data_reader;
class Blas_matrix;


class Self_energy {

  struct Conf {
    unsigned a,k,l,type;
  };

  
  Integral_table& table_;
  SCF_data_reader& phis_;
  
  std::vector<Conf> sats_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> occs_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> virs_[8];

  double V1212(unsigned a, unsigned b, unsigned c, unsigned d);
  
  void sigma_3diagrams(Blas_matrix& sigma);      // An alternative way to compute \Sigma to third order.
  void dynamic_self_energy_3(Blas_matrix& M_ak); // Computes the third-order dynamic self-energy, particle-hole part.
  void linear_eq_selfenergy(Blas_matrix& b);     // Solves a set of non-linear equations for \Sigma4+.

  void rho2sigma(Blas_matrix& rho);              // Converts lower-order density to higher-order self-energy.

  void rho_hole_part_2(Blas_matrix& rho);        
  void rho_holeparticle_2(Blas_matrix& rho);
  void rho_particle_part_2(Blas_matrix& rho);

  void rho_hole_part_3(Blas_matrix& rho);
  void rho_holeparticle_3(Blas_matrix& rho, Blas_matrix& self_energy);
  void rho_particle_part_3(Blas_matrix& rho);
public:
  
  Self_energy(Integral_table& tb, SCF_data_reader& r);

  void density_2(Blas_matrix& rho); // Produces the second-order one-particle density matrix.
  void density_3(Blas_matrix& rho); // Produces the third-order one-particle density matrix.
  void density_3plus(Blas_matrix& rho); // Produces an improved third-order one-particle density matrix,
                                        // making use of the 4+ static self-energy.
  void static_selfenergy_3(Blas_matrix& sigma); // Produces the third order static self-energy.
  void static_selfenergy_4(Blas_matrix& sigma); // Produces the fourth order static self-energy.
  void static_selfenergy_4plus(Blas_matrix& sigma); // Produces the static self-energy using the 4+ scheme.

};



#endif //#ifndef __SELF_ENERGY_HPP__

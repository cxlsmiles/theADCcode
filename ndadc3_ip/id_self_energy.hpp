#ifndef __IDSELF_ENERGY_HPP__
#define __IDSELF_ENERGY_HPP__


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


class idSelf_energy {

  struct Conf {
    unsigned a,k,l,type;
  };



  typedef double (idSelf_energy::*Term)(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  typedef void (idSelf_energy::*Term1)(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);



  double eq20term_1(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_2(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_3(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_4(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_6(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_7(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq20term_12(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);


  double eq_c12term_1(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_2(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_3(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_4(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_6(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_7(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
  double eq_c12term_12(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);


  void eq_c16term_1(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_2(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_3(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_4(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_6(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_7(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_c16term_12(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);




  void eq_30term_1(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_2(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_3(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_4(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_6(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_7(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);
  void eq_30term_12(unsigned, unsigned, unsigned, unsigned, unsigned, Blas_matrix&, double);


  void eq20(Blas_matrix& rho,  Term func, double eri, unsigned a, unsigned b, unsigned c, unsigned d);
  void eq_c12(Blas_matrix& M_ak, Term func, double eri, unsigned a, unsigned b, unsigned c, unsigned d);
  void eq_c16(Blas_matrix& M_ak, Term1 func, double eri, unsigned a, unsigned b, unsigned c, unsigned d);
  void eq_30(Blas_matrix& M_ak, Term1 func, double eri, unsigned a, unsigned b, unsigned c, unsigned d);




  void contrib_4vir(Blas_matrix& rho, Blas_matrix& M_ak, double eri, unsigned a, unsigned c, unsigned d, unsigned b);

  void eq20term(Blas_matrix& rho, double eri1, double eri2,  unsigned a, unsigned c, unsigned d, unsigned b);
  void c12c16term(Blas_matrix& rho, double eri, double eri1, unsigned a, unsigned c, unsigned d, unsigned b);
  void eq30term(Blas_matrix& rho, double eri, unsigned a, unsigned c, unsigned d, unsigned b);



  bool get_next_fourvir(unsigned& vir1, unsigned& vir2, unsigned& vir3, unsigned& vir4, double& vi);


  Integral_table& table_;
  SCF_data_reader& phis_;
  
  std::vector<Conf> sats_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> occs_[8]; // separate the satellite states into symmetries
  std::vector<unsigned> virs_[8];

  double V1212(unsigned a, unsigned b, unsigned c, unsigned d);
  
  void sigma_3diagrams(Blas_matrix& sigma);      // An alternative way to compute \Sigma to third order.
  void dynamic_self_energy_3_res(Blas_matrix& M_ak); // Computes the third-order dynamic self-energy, particle-hole part.


  void linear_eq_selfenergy(Blas_matrix& b);     // Solves a set of non-linear equations for \Sigma4+.

  void rho2sigma(Blas_matrix& rho);              // Converts lower-order density to higher-order self-energy.

  void rho_hole_part_2(Blas_matrix& rho);        
  void rho_holeparticle_2(Blas_matrix& rho);
  void rho_particle_part_2(Blas_matrix& rho);

  void rho_hole_part_3_res(Blas_matrix& rho);
  void rho_holeparticle_3(Blas_matrix& rho, Blas_matrix& self_energy);
  void rho_particle_part_3_res(Blas_matrix& rho);


  // integral driven routines

  



public:
  
  idSelf_energy(Integral_table& tb, SCF_data_reader& r);

  void static_selfenergy_4plus(Blas_matrix& sigma); // Produces the static self-energy using the 4+ scheme.

};



#endif //#ifndef __IDSELF_ENERGY_HPP__

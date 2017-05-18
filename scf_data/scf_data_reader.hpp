#ifndef __SCF_DATA_READER_HPP__
#define __SCF_DATA_READER_HPP__

#include "scf_data_interface.hpp"
#include <vector>




class SCF_data_reader: public Scf_data_interface {
  

  static const unsigned mult_table_[][8];
  unsigned number_irreps_;
  unsigned number_occupied_;
  unsigned number_orbitals_;
  std::vector<int> orbital_symmetry_;
  std::vector<double> orbital_energy_;
 
  void get_info(const char *inp_backend);
  void get_sym();
  void get_epsi();
  void get_occ();

  Scf_data_interface* scf;

public:
  
  SCF_data_reader(const char *inp_backend);
  ~SCF_data_reader();

  
  void get_info (int *nSym, int *nBas, int *nCenters);
  void get_epsi (double *E_hf, double *e, int *n);
  void get_sym (int *s, int *n);
  void get_occ (double *o, int *n);
  void get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint);
  void get_scfvec (double *C, int *n, int *len);
  void get_overlap (double *S, int *n);
  void get_dip (double *x, double *y, double *z, int *n);
  void get_vel (double *x, double *y, double *z, int *n);
  void get_quad(double *xx, double *xy, double *xz,
			     double *yy, double *yz, double *zz, int *n);
  void get_oneel(double *h, double *t, int *n);
  void get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc);
  void get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo);


  double get_integral(int& p, int& q, int& r, int& s);
  inline unsigned number_occupied()  {return number_occupied_;}
  inline unsigned number_irreps() {return number_irreps_;}
  inline unsigned number_orbitals() {return number_orbitals_;}
  inline double energy(unsigned orb)  {return orbital_energy_[orb];}
  //determines the symmetry of a given orbital
  inline unsigned irrep(unsigned orb) {return orbital_symmetry_[orb];}
  static inline unsigned irrep_product(unsigned sym1, unsigned sym2) {return mult_table_[sym1][sym2];}
 
};

#endif //#ifndef __SCF_DATA_READER_HPP__

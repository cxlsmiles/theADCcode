#ifndef __SCF_DATA_INTERFACE_HPP__
#define __SCF_DATA_INTERFACE_HPP__


//see phis/SPEC

class Scf_data_interface {
public:
 
  virtual void get_info (int *nSym, int *nBas, int *nCenters) = 0;
  virtual void get_epsi (double *E_hf, double *e, int *n) = 0;
  virtual void get_sym (int *s, int *n) = 0;
  virtual void get_occ (double *o, int *n) = 0;
  virtual void get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint) = 0;
  virtual void get_scfvec (double *C, int *n, int *len) = 0;
  virtual void get_overlap (double *S, int *n) = 0;
  //Returns the components of the dipole matrix (length form) in much the same
  //way as the overlap integrals are returned.
  virtual void get_dip (double *x, double *y, double *z, int *n) = 0;

  //Returns the components of the dipole matrix (velocity form) just like
  //phis_get_dip() above.
  virtual void get_vel (double *x, double *y, double *z, int *n) = 0;

  //Returns the components of the quadrupole moment.
  virtual void get_quad(double *xx, double *xy, double *xz,
			     double *yy, double *yz, double *zz, int *n) = 0;

  //Returns the square one-electron matrices in a form similar to the overlap
  //integrals above. The one-electron Hamiltonian is returned in the storage
  //area indicated by the first argument. The kinetic integrals are stored in
  //the area passed as the second argument.
  virtual void get_oneel(double *h, double *t, int *n) = 0;

  //uknown
  virtual void get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc) = 0;
  // computes MO cap, implemented only for GUS
  virtual void get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo) = 0;

};

#endif //#ifndef __SCF_DATA_INTERFACE_HPP__

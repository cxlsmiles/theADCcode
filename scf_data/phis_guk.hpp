#ifdef GUK
#ifndef __PHIS_GUK_HPP__
#define __PHIS_GUK_HPP__

#include "scf_data_interface.hpp"

//see phis/SPEC

class Phis_guk: public Scf_data_interface {
public:
  Phis_guk(const char *inp_backend);


  virtual void get_info (int *nSym, int *nBas, int *nCenters);
  virtual void get_epsi (double *E_hf, double *e, int *n);
  virtual void get_sym (int *s, int *n);
  virtual void get_occ (double *o, int *n);
  virtual void get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint);
  virtual void get_scfvec (double *C, int *n, int *len);
  virtual void get_overlap (double *S, int *n);
  virtual void get_dip (double *x, double *y, double *z, int *n);
  virtual void get_vel (double *x, double *y, double *z, int *n);
  virtual void get_quad(double *xx, double *xy, double *xz,
			double *yy, double *yz, double *zz, int *n);
  virtual void get_oneel(double *h, double *t, int *n);
  virtual void get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc);
  virtual void get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo){}
};

#endif //#infdef __PHIS_GUK_HPP__
#endif

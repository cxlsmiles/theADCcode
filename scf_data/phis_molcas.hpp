#ifdef MOLC
#ifndef __PHIS_MOLCAS_HPP__
#define __PHIS_MOLCAS_HPP__

#include "scf_data_interface.hpp"

//see phis/SPEC

#include <vector>

class Phis_molcas: public Scf_data_interface {


  void sort(); 
  std::vector<unsigned> ind;
  std::vector<unsigned> inv_ind;

  template <class T>
  void reorder(T* mat, int dim1, int dim2)
  {
    
    if ((dim1 < 0) || (dim2 < 0)) return;
    
    T* temp_mat = new  T[dim1*dim2];
    
    if (dim2 == 1) dim2 = 0;
    
    
    for (int i = 0; i < dim1; i++) {
      int j = 0;
      do {
	temp_mat[i +j * dim2] = mat[i + j*dim2];
	j++;
      } while (j < dim2);
    }
    
    for (int i = 0; i < dim1; i++) {
      int j = 0;
      do {
	mat[i+j*dim2] = temp_mat[ind[i] + ind[j]*dim2] ;
	j++;
      } while (j < dim2);
    }
    
    
    delete [] temp_mat;
    
  }


    
public:
  Phis_molcas(const char *inp_backend);


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
  virtual void get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo);
  void molcas2capinput();
};

#endif //#infdef __PHIS_MOLCAS_HPP__
#endif

#ifdef GUK
#include "phis_guk.hpp"
#include <string>


extern "C" {
#include "phis_guk.h"
}

#include <iostream>
using namespace std;



Phis_guk::Phis_guk(const char *inp_backend)
{
  string backend_name(inp_backend);
  int backend;
  
  /* choose backend for PHIS */
  if(!backend_name.compare("molcas"))
    backend=MOLCAS | (SYM_BLOCKED << 16);
  else if(!backend_name.compare("guk"))
    backend=GUK;
  else 
    throw string("SCF_data_reader::get_info(): Uknown Backend \
Available: guk, molcas.\n");

  int cap;
  /* initialize PHIS */
  cap=guk_phis_init(&backend);
  
  /* check if backend has the needed functions */
  if(~(cap | ~(HAVE_INFO | HAVE_LOA | HAVE_EPSI | HAVE_SYM | HAVE_OCC))) {

    throw string("SCF_data_reader::get_info: Backend incomplete. \
Terminating.\n");
  } 
}

void Phis_guk::get_info (int *nSym, int *nBas, int *nCenters)
{
  guk_phis_get_info(nSym, nBas, nCenters);
}
void Phis_guk::get_epsi (double *E_hf, double *e, int *n)
{
  guk_phis_get_epsi (E_hf, e, n);
}
void Phis_guk::get_sym (int *s, int *n)
{
  guk_phis_get_sym (s, n);
}
void Phis_guk::get_occ (double *o, int *n)
{
  guk_phis_get_occ (o, n);
}
void Phis_guk::get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint)
{
  guk_phis_get_next_Vpqrs(p, q, r, s, vint);
}
void Phis_guk::get_scfvec (double *C, int *n, int *len)
{
  guk_phis_get_scfvec (C, n, len);
}
void Phis_guk::get_overlap (double *S, int *n)
{
  guk_phis_get_overlap (S, n);
}

void Phis_guk::get_dip (double *x, double *y, double *z, int *n)
{
  guk_phis_get_dip (x, y, z, n);
}
void Phis_guk::get_vel (double *x, double *y, double *z, int *n)
{
  guk_phis_get_vel (x, y, z, n);
}
void Phis_guk::get_quad(double *xx, double *xy, double *xz,
			     double *yy, double *yz, double *zz, int *n)
{
  guk_phis_get_quad(xx, xy, xz, yy, yz, zz, n);
}
void Phis_guk::get_oneel(double *h, double *t, int *n)
{
  guk_phis_get_oneel(h, t, n);
}

void Phis_guk::get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc)
{
  guk_phis_get_geometry(n, geometry, Z_nuc, E_nuc);
}


#endif

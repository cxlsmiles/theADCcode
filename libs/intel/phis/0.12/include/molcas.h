#ifndef MOLCAS_H
#define MOLCAS_H

#include "phis.h"

/* Wrappers for FORTRAN subroutine calls */
#define mc_setup F77_FUNC_(mc_setup,MC_SETUP)
#define mc_query_vel F77_FUNC_(mc_query_vel,MC_QUERY_VEL)
#define mc_stop F77_FUNC_(mc_stop,MC_STOP)

#define phis_mc_info F77_FUNC_(phis_mc_info,PHIS_MC_INFO)
#define phis_mc_epsi F77_FUNC_(phis_mc_epsi,PHIS_MC_EPSI)
#define phis_mc_sym F77_FUNC_(phis_mc_sym,PHIS_MC_SYM)
#define phis_mc_occ F77_FUNC_(phis_mc_occ,PHIS_MC_OCC)
#define phis_mc_nextint F77_FUNC_(phis_mc_nextint,PHIS_MC_NEXTINT)
#define phis_mc_scfvec F77_FUNC_(phis_mc_scfvec,PHIS_MC_SCFVEC)
#define phis_mc_overlap F77_FUNC_(phis_mc_overlap,PHIS_MC_OVERLAP)
#define phis_mc_geometry F77_FUNC_(phis_mc_geometry,PHIS_MC_GEOMETRY)
#define phis_mc_ao F77_FUNC_(phis_mc_ao,PHIS_MC_AO)
#define phis_mc_active F77_FUNC_(phis_mc_active,PHIS_MC_ACTIVE)
#define phis_mc_dip F77_FUNC_(phis_mc_dip,PHIS_MC_DIP)
#define phis_mc_vel F77_FUNC_(phis_mc_vel,PHIS_MC_VEL)
#define phis_mc_quad F77_FUNC_(phis_mc_quad,PHIS_MC_QUAD)
#define phis_mc_oneel F77_FUNC_(phis_mc_oneel,PHIS_MC_ONEEL)

/* These are local helper functions to be called from m_init */
void mc_setup(int *flags, double *work, int *len);
int mc_query_vel(void);
void mc_stop(void);

/* These are the API functions */
void phis_mc_info(int *nSym, int *nBas, int *nCenters);
void phis_mc_epsi(double *E_hf, double *epsi, int *n);
void phis_mc_sym(int *sym, int *n);
void phis_mc_occ(double *occ, int *n);
void phis_mc_nextint(int *i, int *j, int *k, int *l, double *Vpqrs);
void phis_mc_scfvec(double *C, int *n, int *len);
void phis_mc_overlap(double *S, int *n);
void phis_mc_geometry(int *n, double *geo, double *Znuc, double *Enuc);
void phis_mc_ao(int *polynomial, int *nmb_cc, double *cc, double *alpha, 
		int *center, int *nmb_ao, int *max_nmb_cc);
void phis_mc_active(int *list, int *n); 
void phis_mc_dip(double *x, double *y, double *z, int *dim);
void phis_mc_vel(double *x, double *y, double *z, int *dim);
void phis_mc_quad(double *xx, double *xy, double *xz,
		  double *yy, double *yz, double *zz, int *dim);
void phis_mc_oneel(double *h, double *t, int *n);

#endif

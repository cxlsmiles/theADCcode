/* include/phis.h.  Generated from phis.h.in by configure.  */
#ifndef PHIS_H
#define PHIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Implemented packages */
#define MAX_HOOKS 10
#define MOLCAS 0
#define GUK    1
#define STUB   2
#define MAXLENGTH 4096

/* Available capabilities */
#define HAVE_INFO       0x0001
#define HAVE_EPSI       0x0002
#define HAVE_SYM        0X0004
#define HAVE_OCC        0x0008
#define HAVE_INT        0x0010
#define HAVE_NEXTINT    0x0020
#define HAVE_SCFVEC     0x0040
#define HAVE_OVERLAP    0x0080
#define HAVE_GEOMETRY   0x0100
#define HAVE_AO         0x0200
#define HAVE_LOA        0x0400
#define HAVE_DIP        0x0800
#define HAVE_VEL        0x1000
#define HAVE_QUAD	0x2000
#define HAVE_ONEEL	0x4000

/* Supported flags */
#define FLAG_RSVD1      0x0001
#define SYM_BLOCKED     0x0002

/*
 * global variables and useful macros
 */

/* export the multiplication table */
extern int F77_FUNC(multab,MULTAB)[64];

/* use this macro to access multab in C */
#define MULTAB(a,b) (F77_FUNC(multab,MULTAB)[(a-1)*8+b-1])


/*
 * generic interfaces; these are called from applications.
 * C function calls are actually macros to hide FORTRAN calling conventions
 * from the application.
 */

#define phis_init F77_FUNC_(phis_init,PHIS_INIT)
#define phis_get_info F77_FUNC_(phis_get_info,PHIS_GET_INFO)
#define phis_get_epsi F77_FUNC_(phis_get_epsi,PHIS_GET_EPSI)
#define phis_get_sym F77_FUNC_(phis_get_sym,PHIS_GET_SYM)
#define phis_get_occ F77_FUNC_(phis_get_occ,PHIS_GET_OCC)
#define phis_get_next_Vpqrs F77_FUNC_(phis_get_next_vpqrs,PHIS_GET_NEXT_VPQRS)
#define phis_get_scfvec F77_FUNC_(phis_get_scfvec,PHIS_GET_SCFVEC)
#define phis_get_overlap F77_FUNC_(phis_get_overlap,PHIS_GET_OVERLAP)
#define phis_get_geometry F77_FUNC_(phis_get_geometry,PHIS_GET_GEOMETRY)
#define phis_get_ao F77_FUNC_(phis_get_ao,PHIS_GET_AO)
#define phis_list_active F77_FUNC_(phis_list_active,PHIS_LIST_ACTIVE)
#define phis_get_dip F77_FUNC_(phis_get_dip,PHIS_GET_DIP)
#define phis_get_vel F77_FUNC_(phis_get_vel,PHIS_GET_VEL)
#define phis_get_quad F77_FUNC_(phis_get_quad,PHIS_GET_QUAD)
#define phis_get_oneel F77_FUNC_(phis_get_oneel,PHIS_GET_ONEEL)
#define phis_mc_nextint_close F77_FUNC_(phis_mc_nextint_close,PHIS_MC_NEXTINT_CLOSE)

#define phis_load_Vpqrs F77_FUNC_(phis_load_vpqrs, PHIS_LOAD_VPQRS)
#define Vpqrs F77_FUNC(vpqrs,VPQRS)
#define Vordered F77_FUNC(vordered,VORDERED)
#define phis_version F77_FUNC_(phis_version,PHIS_VERSION)
#define dgemm F77_FUNC(dgemm, DGEMM)

/* Initialize the package */

int phis_init(int *module);

/* Get general information */
void phis_get_info(int *nSym, int *nBas, int *nCenters);

/* Get orbital energies */
void phis_get_epsi(double *Ehf, double *e, int *n);

/* Symmetry info */
void phis_get_sym(int *s, int *n);

/* Occupation numbers */
void phis_get_occ(double *o, int *n);

/* next integral on disk, along with its indices */
void phis_get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint);

/* Get the SCF vectors */
void phis_get_scfvec(double *C, int *n, int *len);

/* Get the overlap matrix */
void phis_get_overlap(double *S, int *n);

/* Get nuclear configuration */
void phis_get_geometry(int *n, double *geo, double *Z_nuc, double *E_nuc);

/* Get AOs */
void phis_get_ao(int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
		 int *nmb_ao,int *max_nmb_cc);

/* Map from MOs to active orbitals */
void phis_list_active(int *list, int *n);

/* Get dipole integrals (length form). */
void phis_get_dip(double *x, double *y, double *z, int *n);

/* Get dipole integrals (velocity form). */
void phis_get_vel(double *x, double *y, double *z, int *n);

/* Get quadrupole integrals */
void phis_get_quad(double *xx, double *xy, double *xz,
                   double *yy, double *yz, double *zz, int *n);

/* Get one-electron integrals */
void phis_get_oneel(double *h, double *t, int *n);

/* Close TRAINT file. Analog of rewind function*/
void phis_mc_nextint_close(void);
/* -------------------- */

/* load integrals and initialize lookup */
void phis_load_Vpqrs(void);

/* two-electron integral for a given set of orbital indices */
double Vpqrs(const int *p, const int *q, const int *r, const int *s);

/* two-electron integral, optimized version for pre-ordered indices */
double Vordered(const int *p, const int *q, const int *r, const int *s);

/* Function to query the version number of this library. */
void phis_version(int *major, int *minor, int *patch);

/* malloc wrapper, used in some C functions. */
void *malloc_we(size_t size);

/* number of PHIS's variables */
int num_projects(void);

/*list of PHIS's variables */
void list_projects(int *nump,char *list);


/* BLAS general purpose matrix multiplication */
void dgemm(char *transa, char*transb, int *m, int*n, int *k, double *alpha,
	   double *A, int *lda, double *B, int *ldb, 
	   double *beta, double *C, int *ldc);

/* 
 * This is filled in by the init function of the backend.
 */
struct entry_points {
	int cap;
	void (*info) (int *nSym, int *nBas, int *nCenters);
	void (*epsi) (double *E_hf, double *epsi, int *n);
	void (*sym) (int *sym, int *n);
	void (*occ) (double *occ, int *n);
	void (*next_Vpqrs) (int *p, int *q, int *r, int *s, double *Vpqrs);
	void (*scfvec) (double *C, int *n, int *len);
	void (*overlap) (double *S, int *n);
	void (*geometry) (int *, double *, double *, double *);
	void (*ao) (int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
		    int *nmb_ao,int *max_nmb_cc);
	void (*list_active) (int *list, int *n);
	void (*dip) (double *x, double *y, double *z, int *dim);
	void (*vel) (double *x, double *y, double *z, int *dim);
	void (*quad) (double *xx, double *xy, double *xz,
                      double *yy, double *yz, double *zz, int *n);
	void (*oneel) (double *h, double *t, int *n);
	void *unusedE;
	void *unusedF;
};

struct entry_points interface;

#endif

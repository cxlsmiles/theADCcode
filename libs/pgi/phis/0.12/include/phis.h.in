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
extern int F77_FUNC(multabguk,MULTABGUK)[64];

/* use this macro to access multab in C */
#define MULTAB(a,b) (F77_FUNC(multabguk,MULTABGUK)[(a-1)*8+b-1])


/*
 * generic interfaces; these are called from applications.
 * C function calls are actually macros to hide FORTRAN calling conventions
 * from the application.
 */

#define guk_phis_init F77_FUNC_(guk_phis_init,GUK_PHIS_INIT)
#define guk_phis_get_info F77_FUNC_(guk_phis_get_info,GUK_PHIS_GET_INFO)
#define guk_phis_get_epsi F77_FUNC_(guk_phis_get_epsi,GUK_PHIS_GET_EPSI)
#define guk_phis_get_sym F77_FUNC_(guk_phis_get_sym,GUK_PHIS_GET_SYM)
#define guk_phis_get_occ F77_FUNC_(guk_phis_get_occ,GUK_PHIS_GET_OCC)
#define guk_phis_get_next_Vpqrs F77_FUNC_(guk_phis_get_next_vpqrs,GUK_PHIS_GET_NEXT_VPQRS)
#define guk_phis_get_scfvec F77_FUNC_(guk_phis_get_scfvec,GUK_PHIS_GET_SCFVEC)
#define guk_phis_get_overlap F77_FUNC_(guk_phis_get_overlap,GUK_PHIS_GET_OVERLAP)
#define guk_phis_get_geometry F77_FUNC_(guk_phis_get_geometry,GUK_PHIS_GET_GEOMETRY)
#define guk_phis_get_ao F77_FUNC_(guk_phis_get_ao,GUK_PHIS_GET_AO)
#define guk_phis_list_active F77_FUNC_(guk_phis_list_active,GUK_PHIS_LIST_ACTIVE)
#define guk_phis_get_dip F77_FUNC_(guk_phis_get_dip,GUK_PHIS_GET_DIP)
#define guk_phis_get_vel F77_FUNC_(guk_phis_get_vel,GUK_PHIS_GET_VEL)
#define guk_phis_get_quad F77_FUNC_(guk_phis_get_quad,GUK_PHIS_GET_QUAD)
#define guk_phis_get_oneel F77_FUNC_(guk_phis_get_oneel,GUK_PHIS_GET_ONEEL)
#define guk_phis_mc_nextint_close F77_FUNC_(guk_phis_mc_nextint_close,GUK_PHIS_MC_NEXTINT_CLOSE)

#define guk_phis_load_Vpqrs F77_FUNC_(guk_phis_load_vpqrs, GUK_PHIS_LOAD_VPQRS)
#define guk_Vpqrs F77_FUNC(guk_vpqrs,GUK_VPQRS)
#define guk_Vordered F77_FUNC(vordered,VORDERED)
#define guk_phis_version F77_FUNC_(guk_phis_version,GUK_PHIS_VERSION)
#define dgemm F77_FUNC(dgemm, DGEMM)

/* Initialize the package */

int guk_phis_init(int *module);

/* Get general information */
void guk_phis_get_info(int *nSym, int *nBas, int *nCenters);

/* Get orbital energies */
void guk_phis_get_epsi(double *Ehf, double *e, int *n);

/* Symmetry info */
void guk_phis_get_sym(int *s, int *n);

/* Occupation numbers */
void guk_phis_get_occ(double *o, int *n);

/* next integral on disk, along with its indices */
void guk_phis_get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint);

/* Get the SCF vectors */
void guk_phis_get_scfvec(double *C, int *n, int *len);

/* Get the overlap matrix */
void guk_phis_get_overlap(double *S, int *n);

/* Get nuclear configuration */
void guk_phis_get_geometry(int *n, double *geo, double *Z_nuc, double *E_nuc);

/* Get AOs */
void guk_phis_get_ao(int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
		 int *nmb_ao,int *max_nmb_cc);

/* Map from MOs to active orbitals */
void guk_phis_list_active(int *list, int *n);

/* Get dipole integrals (length form). */
void guk_phis_get_dip(double *x, double *y, double *z, int *n);

/* Get dipole integrals (velocity form). */
void guk_phis_get_vel(double *x, double *y, double *z, int *n);

/* Get quadrupole integrals */
void guk_phis_get_quad(double *xx, double *xy, double *xz,
                   double *yy, double *yz, double *zz, int *n);

/* Get one-electron integrals */
void guk_phis_get_oneel(double *h, double *t, int *n);

/* Close TRAINT file. Analog of rewind function*/
void phis_mc_nextint_close(void);
/* -------------------- */

/* load integrals and initialize lookup */
void guk_phis_load_Vpqrs(void);

/* two-electron integral for a given set of orbital indices */
double guk_Vpqrs(const int *p, const int *q, const int *r, const int *s);

/* two-electron integral, optimized version for pre-ordered indices */
double guk_Vordered(const int *p, const int *q, const int *r, const int *s);

/* Function to query the version number of this library. */
void guk_phis_version(int *major, int *minor, int *patch);

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

struct entry_points guk_interface;

#endif

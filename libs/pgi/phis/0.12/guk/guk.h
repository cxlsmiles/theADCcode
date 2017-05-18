#ifndef GUK_H
#define GUK_H

#include <stdio.h>
#include <stdlib.h>

#include "./typedefs.h"

#define FALSE 0
#define TRUE  1

enum type2 { OLAP = 0, KIN, CORE, DIPX, DIPY, DIPZ, ENUC };

/* dfile specification defined in guk_init */
extern int guk_mxprim;
extern int guk_maxaqm;
extern int guk_maxat;
extern int guk_mxshel;
extern int guk_maxorb;

/* Make data availale everywhere in the backend */
extern Type_0 Summary;
extern Type_3 *section3;
extern Type_51 *section51;

/*extern char *dumpfile_name;*/
extern int nBas, nAct, nAtoms;
extern int *map_active;

/* guk_get_block passes data via a global array */
extern char Guk_buffer[4096];

/* Hard-coded names of ed3 && ed6*/
static char *dumpfile_name="dfile";
static char *mainfile_name="vfile";

/* These are the API functions */
void phis_guk_info(int *nSym, int *nBas, int *nCenters);
void phis_guk_epsi(double *E_hf, double *epsi, int *n);
void phis_guk_sym(int *sym, int *n);
void phis_guk_occ(double *occ, int *n);
void phis_guk_next_Vpqrs(int *p, int *q, int *r, int *s, double *Vpqrs);
void phis_guk_scfvec(double *C, int *n, int *stride);
void phis_guk_overlap(double *S, int *n);
void phis_guk_ao(int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
		 int *nmb_ao,int *max_nmb_cc);
void phis_guk_loa(int *list, int *n);
void phis_guk_geometry(int *n,double *geo,double *Z_nuc,double *e_nuc);
void phis_guk_dip(double *x, double *y, double *z, int *n);

/* These are helper routines in guk_aux.c */
FILE *guk_open_dump_file( char *filename );
void guk_position_file( FILE *fp, int pos );
int guk_get_block( FILE *fp, int blkno );
void *guk_get_section( FILE *fp, int type );
int guk_read_data( FILE *fp, char *buf, int blkno, int n_words );
int guk_sec_by_type( int type );

/* helper functions in typeinits.c */
Type_1* guk_init_Type_1(Type_1* sec);
Type_3* guk_init_Type_3(Type_3* sec);
Type_51* guk_init_Type_51(Type_51* sec);
Type_1005* guk_init_Type_1005(Type_1005* sec);
char* guk_init_section(char* section, int type);
int guk_sizeof_subsection(FILE *fp);
void guk_getpara(FILE *dump_fp);

/* more helper functions */
void sort_gukscf(double *gukscf,Type_3 *gukvec);
void copy_and_expand(double *M, double *buf, int n);
void guk5_read_type2(double *buf, int buflen, int what);
void guk6_read_type2(double *buf, int buflen, int what);
void (*guk_read_type2)(double *buf, int buflen, int what);

#endif

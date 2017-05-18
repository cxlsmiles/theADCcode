/* $Id: globals.h,v 1.6 2003/12/01 14:54:55 joergb Exp $
 *
 * $Log: globals.h,v $
 * Revision 1.6  2003/12/01 14:54:55  joergb
 * parameters changed (make_symtab does not need the
 * whole inp structure => usable in other programs)
 *
 * Revision 1.5  2003/10/05 08:49:46  joergb
 * splitted into globals.h, adc_macros.h, get_input.h,
 *
 * Revision 1.4  2003/08/19 16:55:07  u85
 * new macro needed for satellites: FOR_ALL_2H1P_AKL(A,B),
 * prototypes for the new functions
 *
 * Revision 1.3  2003/03/03 14:23:54  u85
 * some variables added
 *
 * Revision 1.2  2003/01/10 16:31:06  u85
 * f_matrix added and bug in c11_3 (h.c. part) fixed
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */


#include "get_input.h"
#include "read_scf_data.h"
#include "make_symtab.h"
#include "adc_macros.h"


#ifdef DEF_GLOBALS
#define EXTERN 
#else
#define EXTERN extern
#endif


void get_input(char *defaults,int parse_stdin,struct inp_t *inp);
void init_files(struct inp_t *inp,struct symtab_t *symtab);
//void read_scf_data(struct inp_t *inp, struct scf_t *scf);
void read_scf_data(char *backend, int debug, FILE *data_fp, struct scf_t *scf);
void print_sum_formula(int nAtoms,double *Z_nuc);
void make_symtab(int sym,int affinity_mode,int debug,FILE *data_fp,
                 struct scf_t *scf,struct symtab_t *symtab);
void make_matrix(int rows,int cols,char mode,char fill,double **matrix);

int calc_dim_2h1p(struct inp_t *inp,struct symtab_t *symtab);
int calc_dim_2p1h(struct inp_t *inp,struct symtab_t *symtab);

void build_nd_matrices(struct inp_t *inp,struct scf_t *scf,
		      struct symtab_t *symtab);

void calc_k1(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
	     double *c_matrix,double *f_matrix);
void calc_c11_2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix,double *f_matrix);
void calc_c11_3(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix,double *f_matrix);

void calc_c12_1(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix);
void calc_c12_2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix);

void calc_k2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
	     double *c_diag);
void calc_c22_1_diag(struct inp_t *inp,struct scf_t *scf,
		     struct symtab_t *symtab,double *c_diag);
void calc_c22_1_cols(struct inp_t *inp,struct scf_t *scf,
		     struct symtab_t *symtab,
		     int *m,int *n,int *b,int row,double *c_col0,double *c_col1);
void calc_c22_1_off(struct inp_t *inp,struct scf_t *scf,
		    struct symtab_t *symtab,int dim_2h1p,double *c_22);

void convert_dyson_matrix(struct inp_t *inp,struct symtab_t *symtab,
			  double *hh_matrix);
void convert_sigma_matrix(struct inp_t *inp,struct symtab_t *symtab,
			  double *hh_matrix);

void write_matrices(struct inp_t *inp,struct scf_t *scf,
		    struct symtab_t *symtab,
		    double *c11,double *c12,double *c22_diag,double *f);

/* These functions supply Joerg's code with the
   two electron integrals and symmetry information. */
double vpqrs_(int*, int*, int*, int*);
int    Multab_(int, int);

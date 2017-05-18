/* $Id: globals.h,v 1.0 date time imke Exp$
 *
 * $Log: globals.h,v $
 * Revision 1.1 time date imke  14/06/04
 * new function  read_all_cap_elements 
 */

#include "ndadc3_ip/globals.h" 
#include "my_adc_macros.h"

#ifdef DEF_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif


/*void get_input_cap(char *defaults, int parse_stdin, struct inp_t *inp);
void cap_init_files(struct inp_t *inp,struct symtab_t *symtab); 
void read_scf_data(struct inp_t *inp,struct scf_t *scf);*/
/* void make_symtab(int sym,int affinity_mode,int debug,FILE *data_fp,  */
/*                  struct scf_t *scf, struct symtab_t *symtab); */
/* void make_matrix(int rows,int cols,char mode,char fill,double **matrix); */

/* int calc_dim_2h1p(struct inp_t *inp,struct symtab_t *symtab); */
/* int calc_dim_2p1h(struct inp_t *inp,struct symtab_t *symtab); */

/* void read_cap_elements(struct inp_t *inp,struct symtab_t *symtab,  */
/*                        double *cap_elements); */

/* void read_all_cap_elements(struct inp_t *inp, struct symtab_t *symtab,  */
/*                        double *cap_elements); */

/* void cap_calc_psi0_rho_psi0(struct inp_t *inp,struct scf_t *scf,  */
/*                         struct symtab_t *symtab,double *rho); */

void cap_calc_d_null(struct scf_t *scf,double *rho_0,double *cap_elements,
                     double *d_null, int debug);

void calc_d_null_null(struct symtab_t *symtab,double *cap_elements,
                              double *d_null_null, int debug);

/* void build_d_matrices(struct inp_t *inp,struct scf_t *scf,  */
/*                       struct symtab_t *symtab); */

void cap_calc_d11(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
                  double *d_null,double *d11_matrix,double *rho_0, 
                  double *cap_elements);

void cap_calc_d12(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
                  double *d12_matrix,double *cap_elements);

void cap_calc_d22_diag(struct inp_t *inp,struct scf_t *scf, 
                       struct symtab_t *symtab,double *d_null_null, 
                       double *d22_diag,double* cap_elements);

void cap_calc_d22_off_cols(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
                           int *m,int *n,int *b,int row,double *d_col0,double *d_col1, 
                           double *cap_elements);

/* void write_d_matrices(struct inp_t *inp,struct scf_t *scf, */
/*                      struct symtab_t *symtab, */
/*                      double *d11,double *d12,double *d22_diag,double *cap_elements); */

/* void print_sum_formula(int nAtoms,double *Z_nuc); */

/*double Vpqrs_cap(int* a, int* b, int* c, int* d);
int Multab_cap(int a, int b);
*/

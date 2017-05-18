/* $Id: get_input.h,v 1.3 2004/01/29 18:47:58 joergb Exp $
 *
 * $Log: get_input.h,v $
 * Revision 1.3  2004/01/29 18:47:58  joergb
 * added F-matrix (identical to prev. Rel.)
 *
 * Revision 1.1  2003/10/05 08:49:46  joergb
 * splitted into globals.h, adc_macros.h, get_input.h,
 *
 */


/*#include "filetools.h"*/


/* inp_t:
 * structure for all input variables
 */
struct inp_t {
  char   *backend;       /* PHIS backend */

  /* FIXME: will vanish */
  char   *store_c;       /* mode to store f-matrix */
  char   *store_f;       /* mode to store f-matrix */
  char   *sigma_in;      /* mode ~adc1p_sigma: sigma output file */
  char   *sigma_out;     /* mode ~adc1p_sigma: sigma output file */

  int    debug;          /* debug level */
/*   FILE   *data_fp;       /\* file pointer for data output *\/ */
  int    sym;            /* problem symmetry */
  int    affinity_mode;  /* switch to ND affinity calculation (experimental) */

  /* inputfiles */
  struct fileinfo_t *c11_file;  /* main/main file */
  struct fileinfo_t *c12_file;  /* main/sat  file */
  struct fileinfo_t *c22_file;  /*  sat/sat  file */
  struct fileinfo_t *cout_file; /* output file */
  struct fileinfo_t *fout_file; /* output file */
  int *c11_file_type;           /* filezypes for these files */
  int *c12_file_type;
  int *c22_file_type;
  int *cout_file_type;
  int *fout_file_type;

  /* orders of the several blocks */
  int calc_c11_orders[4];
  int calc_c12_orders[4];
  int calc_c22_orders[4];
};



void get_input(char *defaults,int parse_stdin,struct inp_t *inp);



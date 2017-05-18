#include "my_globals.h"
/* #include "phis.h" */
#include "my_adc_macros.h"


static char *cvsid="$Id: calc_d_null_null.c,v 1.0 date time imke Exp$";


/* calc_d_null_null.c:
 *
 * calculates the zeroth order operator moment D0(0) needed in ISR formulae
 *
 * D0=sum over all d_k of occupied orbital indeces k
 *               
 *
 * 'INPUT': cap_elements (matrix)
 * 'output: d_null_null number! (double)
 */

void calc_d_null_null(struct symtab_t *symtab,double *cap_elements,double *d_null_null, int debug)
{

  int k_sym, *k;


  for(k_sym=1; k_sym<=symtab->nSym; k_sym++){
    FOR_OCC_IN_SYM(k,k_sym){
      *d_null_null+=CAP_PTR(k,k);
    }
  }

/*   printf("finished calculating d_null_null= %16.14f \n", *d_null_null); */
/*   printf("leaving routine calc_d_null_null ....\n"); */
}



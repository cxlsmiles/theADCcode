/* #include "phis.h" */
#include "my_globals.h"
/* #include "my_adc_macros.h" */

static char *cvsid="$Id: cap_calc_d22_diag.c,v 1.1 22.04.04 18:24:01 imke Exp$";

/* cap_calc_d22_diag:
 *
 * calculates the diagonal Elements of the 2h1p/2h1p-block of the ISR
 * of a one-particle operator (e.g. a CAP)
 * using formulae by Trofimov, Schirmer
 *
 *
 *
 * 'OUTPUT': matrix d22_diag
 * 'INPUT': d_null_null, number ! should be calc. globally, also used for d11 
 * 'INPUT': cap_elements, matrix containing the one-particel-op. matrix els.
 *           (MO/CAP/MO)
 *
 * $Log cap_calc_d22_diag.c,v
 * 1.2 19.01.05 corrected version: commented 1d,1f which are not to appear here
 * 1.1 22.04 imke
 *     correction of sign mistake in k==l part of loop 1c/1e += -> -=
 *
 * 1.0 initial version
 */

/*-----------------------------------------------------------------*/
/* NOMENCLATURE, this time only one set of indeces is needed       */
/* because of delta_aa' delta_ii' delta_jj' (delta_aa' delta_kk'   */
/* delta_ll' ) because we only treat the diagonal here!!!          */

void cap_calc_d22_diag(struct inp_t *inp, struct scf_t *scf, 
                       struct symtab_t *symtab, double *d_null_null, 
                       double *d22_diag, double* cap_elements)
{
  int a_sym, k_sym, l_sym, kI_sym;            /* used for macro for_all... */
  int *a, *k, *l;                             /* used for macro for_all... */

  double sum_1, sum_2;

  double *d_ptr=d22_diag;                      /* goes through the array */

  int i;
  int dim_2h1p;

/*   printf("calculating the 0th order diagonal elements of the ISR of\n"); */
/*   printf("a one-particle operator (e.g. CAP ... )\n"); */
/*   printf("debug: d_null_null= %f \n", *d_null_null); */


  FOR_ALL_2H1P_AKL({
    /* akl part of 2h1p loop (k!=l) */

    /* initialize helper sums */
    sum_1 = sum_2 = 0.0;

    /* (1a) */
    sum_1 += *d_null_null;
    sum_2 += *d_null_null;

    /* (1b) */
    sum_1 += CAP_PTR(a,a);    /*d_element(a,a,cap_elements);*/
    sum_2 += CAP_PTR(a,a);    /*d_element(a,a,cap_elements);*/


    /* (1c) */
    sum_1 -= CAP_PTR(l,l);    /*d_element(l,l,cap_elements);*/
    sum_2 -= CAP_PTR(l,l);    /*d_element(l,l,cap_elements);*/


    /* (1e) */
    sum_1 -= CAP_PTR(k,k);    /*d_element(k,k,cap_elements);*/
    sum_2 -= CAP_PTR(k,k);    /*d_element(k,k,cap_elements);*/


    /* minus from spin-free representation */
    /* the kronecker deltas of these guys CANNOT be fufilled in the diagonal part!*/
    /* (1d) */
    /*sum_1 += -CAP_PTR(k,l); */  /*d_element(k,l,cap_elements);*/
    /*sum_2 +=  CAP_PTR(k,l); */  /*d_element(k,l,cap_elements);*/

    /* minus from spin-free representation */
    /* the kronecker deltas of these guys CANNOT be fufilled in the diagonal part!*/
    /* (1f) */
    /*sum_1 += -CAP_PTR(l,k); */  /*d_element(l,k,cap_elements);*/
    /*sum_2 +=  CAP_PTR(l,k); */  /*d_element(l,k,cap_elements);*/

    /* two matrixentries, go to next entry for each of them  ... */
    *(d_ptr++) += sum_1;
    *(d_ptr++) += sum_2;
  },{
    /* akl part of 2h1p loop (k==l) */

    /* initialize helper sum */
    sum_1 = 0.0;
  
    /* (1a) */
    sum_1 += *d_null_null;

    /* (1b) */
    sum_1 += CAP_PTR(a,a);    /*d_element(a,a,cap_elements);*/

    /* (1c) */
    sum_1 -= CAP_PTR(l,l);    /*d_element(l,l,cap_elements);*/

    /* (1e) */
    sum_1 -= CAP_PTR(k,k);    /*d_element(k,k,cap_elements);*/

    /* (1d), (1f) are 0 due to spinfree representation  */

    /* matrix entry, go to next matrix element ... */
    *(d_ptr++) += sum_1;

  });

/*   if(inp->debug >3){ */
/*     dim_2h1p=calc_dim_2h1p(inp,symtab); */
/*     printf("\n printing diagonal elements of 2h1p/2h1p-block (indeces starting from 1 ...\n"); */
/*     for(i=0; i<dim_2h1p; i++){ */
/*       printf("Index: %i, Element: %16.14f\n", i+1, d22_diag[i]); */
/*     } */
/*     printf("\n\n"); */
/*   } */

}




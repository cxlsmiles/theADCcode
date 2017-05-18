#include "my_globals.h"
/* #include "phis.h" */
#include "my_adc_macros.h"

static char *cvsid="$Id: cap_calc_d_null.c,v 1.0 date time imke Exp$";

/* cap_calc_d_null.c:
 *
 * calculates the second order operator moment D0(2) needed in ISR formulae
 *
 * D0=Tr(d*rho) :d MO representation of one-particle operator
 *              :rho ground state one particel matrix
 *
 * such that d_null(double) = sum_k sum_i d_ki rho_ik 
 *
 * 'INPUT': rho_0, cap_elements (matrices)
 * 'output: d_null number! (double)
 */

void cap_calc_d_null(struct scf_t *scf, double *rho_0, double *cap_elements,
                     double *d_null, int debug)
{
  int i, k;


  for(k=1; k<=scf->nBas; k++){
    for(i=1; i<=scf->nBas; i++){
      *d_null += CAP_INT(k,i) * RHO_0_INT(i,k);    /*element(k,i,cap_elements)*element(i,k,rho_0);*/
      if(debug >2){
	printf(" k: %i i: %i cap_matrix_element: %16.14f rho_null: %16.14f\n", k, i, CAP_INT(k,i), RHO_0_INT(i,k));
      }
    }
  }
   
/*   printf("calculating d_null ....\n");   */
/*   printf("d_null= %16.14f\n", *d_null);  */
/*   printf("leaving routine calc_d_null ...\n\n\n"); */
 
}

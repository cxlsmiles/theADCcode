
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_c22_1_dia.c,v 1.1 2003/08/19 16:04:19 u85 Exp $";

/* calc_c22_1_diag
 *
 * calculates the submatrix K2 (the 0th order of the 2h1p/2h1p-block) 
 *
 * $Log: calc_c22_1_dia.c,v $
 * Revision 1.1  2003/08/19 16:04:19  u85
 * Initial import
 *
 *
 */
void calc_c22_1_diag(struct inp_t *inp,struct scf_t *scf,
		     struct symtab_t *symtab,double *c_diag)
{
  /*  int a_sym,k_sym,l_sym,ki_sym;
      int *a,*k,*l;*/
  double *c_ptr=c_diag;
  double sgn;

  int k_sym,l_sym,a_sym,kI_sym,*k,*l,*a;

/*   printf("calculating 1st order diagonal elements of c22 (sat/sat-block)\n"); */

  /* this term changes sign in affinity mode */
  sgn = inp->affinity_mode ? -1 : 1 ;
  
  FOR_ALL_2H1P_AKL({
    *(c_ptr++) += sgn * 
      ( V1212(k,l,k,l) + V1212(k,l,l,k) 
	- V1212(a,l,a,l) + 0.5 * V1212(a,l,l,a)
	- V1212(a,k,a,k) + 0.5 * V1212(a,k,k,a));
    
    *(c_ptr++) += sgn *
      ( V1212(k,l,k,l) - V1212(k,l,l,k) 
	- V1212(a,l,a,l) + 1.5 * V1212(a,l,l,a)
	- V1212(a,k,a,k) + 1.5 * V1212(a,k,k,a));
  },{
    *(c_ptr++) += sgn *
      ( V1212(k,l,k,l) 
	- V1212(a,l,a,l) + 0.5 * V1212(a,l,l,a)
	- V1212(a,k,a,k) + 0.5 * V1212(a,l,l,a));
  });


#ifdef OLD
  for(k_sym=1;k_sym<=symtab->nSym;k_sym++) {
    ki_sym=MULTAB(inp->sym,k_sym);
    for(l_sym=1;l_sym<k_sym;l_sym++) {
      a_sym=MULTAB(ki_sym,l_sym);

      FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k) {
	FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {
	  if(k_sym==l_sym && *l>=*k) break;
	  FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {

	    c_old=*c_ptr;
	    *(c_ptr++) += sgn * 
	      ( V1212(k,l,k,l) + V1212(k,l,l,k) 
		- V1212(a,l,a,l) + 0.5 * V1212(a,l,l,a)
		- V1212(a,k,a,k) + 0.5 * V1212(a,k,k,a));

	    c_old=*c_ptr;
	    *(c_ptr++) += sgn *
	      ( V1212(k,l,k,l) - V1212(k,l,l,k) 
		- V1212(a,l,a,l) + 1.5 * V1212(a,l,l,a)
		- V1212(a,k,a,k) + 1.5 * V1212(a,k,k,a));
          }
        }
      }
    }

    /* here k and l have the same symmetry => a has input sym */
    l_sym=k_sym;
    a_sym=inp->sym;

    FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k) {
      for(l=symtab->occ[l_sym];l<k;l++) {
	epsi_l=scf->epsi[*l];
	FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {
	  
	  c_old=*c_ptr;
	  *(c_ptr++) += sgn * 
	    ( V1212(k,l,k,l) + V1212(k,l,l,k) 
	      - V1212(a,l,a,l) + 0.5 * V1212(a,l,l,a)
	      - V1212(a,k,a,k) + 0.5 * V1212(a,k,k,a));
	  
	  c_old=*c_ptr;
	  *(c_ptr++) += sgn *
	    ( V1212(k,l,k,l) - V1212(k,l,l,k) 
	      - V1212(a,l,a,l) + 1.5 * V1212(a,l,l,a)
	      - V1212(a,k,a,k) + 1.5 * V1212(a,k,k,a));
	}
      }

      /* at least the akk elements */

      l=k;
      FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {
	c_old=*c_ptr;
	*(c_ptr++) += sgn *
	  ( V1212(k,l,k,l) //+ V1212(k,l,l,k) 
	    - V1212(a,l,a,l) + 0.5 * V1212(a,l,l,a)
	    - V1212(a,k,a,k) + 0.5 * V1212(a,l,l,a));
      }
    }
  }

#endif
  /* FIXME: debug output */
}


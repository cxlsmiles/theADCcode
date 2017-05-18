
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_k2.c,v 1.1 2003/08/19 16:04:19 u85 Exp $";

/* calc_k1:
 *
 * calculates the submatrix K2 (the 0th order of the 2h1p/2h1p-block) 
 *
 * $Log: calc_k2.c,v $
 * Revision 1.1  2003/08/19 16:04:19  u85
 * Initial import
 *
 *
 */
void calc_k2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
	     double *c_diag)
{
  /*  int a_sym,k_sym,l_sym,ki_sym;
      int *a,*k,*l;*/
  double epsi_a,epsi_k,epsi_l;
  double *c_ptr=c_diag;

  int k_sym,l_sym,a_sym,kI_sym,*k,*l,*a;

/*   printf("calculating 0th order diagonal elements of c22 (sat/sat-block)\n"); */

  FOR_ALL_2H1P_AKL({
    epsi_a=scf->epsi[*a];
    epsi_k=scf->epsi[*k];
    epsi_l=scf->epsi[*l];

    *(c_ptr++) = +epsi_a - epsi_k - epsi_l;
    *(c_ptr++) = +epsi_a - epsi_k - epsi_l;
  },{
    epsi_a=scf->epsi[*a];
    epsi_k=scf->epsi[*k];

    *(c_ptr++) = epsi_a - 2. * epsi_k;
  });

#ifdef OLD
  for(k_sym=1;k_sym<=symtab->nSym;k_sym++) {
    ki_sym=MULTAB(inp->sym,k_sym);
    for(l_sym=1;l_sym<k_sym;l_sym++) {
      a_sym=MULTAB(ki_sym,l_sym);

      FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k) {
	FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {
	  if(l_sym==k_sym && *l>=*k) break;
	  FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {

	    *(c_ptr++) = +epsi_a - epsi_k - epsi_l;
	    *(c_ptr++) = +epsi_a - epsi_k - epsi_l;
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
	  
	  *(c_ptr++) = epsi_a - epsi_k - epsi_l;
	  *(c_ptr++) = epsi_a - epsi_k - epsi_l;
	}
      }

      /* at least the akk elements */

      FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {
	*(c_ptr++) = epsi_a - 2. * epsi_k;
      }
    }
  }

  /* FIXME: debug output */

#endif
}


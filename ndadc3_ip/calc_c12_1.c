
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_c12_1.c,v 1.1 2003/08/19 16:04:19 u85 Exp $";

/* calc_c12_1:
 *
 * calculates the submatrix C12_1 (the 1st order of the h/2h1p-block) 
 *
 * Spin-free implementation of formula (A10) [JCP 109, 4734 ('98)], see
 * (A3.3),(A3.6) in the PhD thesis of O. Walter. (A3.3b) is negated here.
 *
 * For the ordering of the coefficients in the 2h1p-block see globals.h
 *
 * $Log: calc_c12_1.c,v $
 * Revision 1.1  2003/08/19 16:04:19  u85
 * Initial import
 *
 *
 */
void calc_c12_1(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix)
{
  int a_sym,k_sym,l_sym,kI_sym;
  int *j,*a,*k,*l;
  double *c_ptr=c_matrix;
  
/*   printf("calculating 1st order of c12 (main/sat-block)\n"); */

  FOR_ALL_2H1P_AKL({
    /* akl part of 2h1p loop (k!=l) */
    FOR_OCC_IN_SYM(j,inp->sym) 
      *(c_ptr++) += SQRT_1_2*(V1212(j,a,k,l)+V1212(j,a,l,k));
    FOR_OCC_IN_SYM(j,inp->sym) 
      *(c_ptr++) -= SQRT_3_2*(V1212(j,a,k,l)-V1212(j,a,l,k)); 
  },{ 
    /* akk part of 2h1p loop */
    FOR_OCC_IN_SYM(j,inp->sym) 
      *(c_ptr++) += V1212(j,a,k,k);
  });

  if(inp->debug<10) return;




  /* debug output follows */

  c_ptr=c_matrix;
  
/*   fprintf(inp->data_fp, */
/* 	  "#<c12-1>\n\n" */
/* 	  "#  a   k   l spin|  j  |  C12(j,akl)\n" */
/* 	  "#----------------+-----+--------------\n"); */
/*   FOR_ALL_2H1P_AKL({ */
/*     FOR_OCC_IN_SYM(j,inp->sym)  */
/*       fprintf(inp->data_fp,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,1,*j,*(c_ptr++)); */
/*     FOR_OCC_IN_SYM(j,inp->sym)  */
/*       fprintf(inp->data_fp,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,2,*j,*(c_ptr++)); */
/*   },{ */
/*     FOR_OCC_IN_SYM(j,inp->sym)  */
/*       fprintf(inp->data_fp,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,1,*j,*(c_ptr++)); */
/*   }); */
/*   fprintf(inp->data_fp,"\n#</c12-1>\n\n"); */
}

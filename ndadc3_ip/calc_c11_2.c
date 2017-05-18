 
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_c11_2.c,v 1.7 2004/01/29 18:48:49 joergb Exp $";

/* calc_c11_2:
 *
 * calculates the submatrix C11_2 (the 2nd order of the hh-block) and the 
 * hh-block of f1_2
 *
 * $Log: calc_c11_2.c,v $
 * Revision 1.7  2004/01/29 18:48:49  joergb
 * corrected formula for hermitian F-matrix (identical to prev. Rel.)
 *
 * Revision 1.5  2003/08/19 16:50:52  u85
 * small changes in output
 *
 * Revision 1.4  2003/06/12 13:42:35  u85
 * debug output changed
 *
 * Revision 1.3  2003/03/03 14:00:57  u85
 * Just prints matrix in debugmode now
 *
 * Revision 1.2  2003/01/10 16:31:04  u85
 * f_matrix added and bug in c11_3 (h.c. part) fixed
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */
void calc_c11_2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix,double *f_matrix)
{
  int a_sym,b_sym,l_sym,ai_sym;
  int *i,*j,*l,*a,*b;
  double epsi_i,epsi_j,epsi_a,epsi_b,epsi_l,epsi_mean,epsi_abl;
  double v_abil,v_abli,v_abjl,v_ablj;
  double vv,c_ij,f_ij;
  int m_ind;
  
  m_ind=0;

/*   printf("calculating 2nd order of c11 (main/main-block)\n"); */

/*   if(inp->debug>9) { */
/*     fprintf(inp->data_fp,"#<c11-2>\n\n" */
/* 	    "#  i   j | C11(i,j)\n" */
/* 	    "#--------+-----------\n"); */
/*   } */

  for(i=symtab->occ[inp->sym];*i!=-1;i++) {
    epsi_i=scf->epsi[*i];

    for(j=symtab->occ[inp->sym];j<=i;j++) {

      epsi_j=scf->epsi[*j];
      epsi_mean=0.5*(epsi_i+epsi_j);

      c_ij=0.;
      f_ij=0.;

      for(a_sym=1;a_sym<=symtab->nSym;a_sym++) {
	ai_sym=MULTAB(inp->sym,a_sym);
	for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	  l_sym=MULTAB(ai_sym,b_sym);
	  
	  FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {
	    FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {
	      FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {

		v_abil=V1212(a,b,i,l);
		v_abli=V1212(a,b,l,i);
		v_abjl=V1212(a,b,j,l);
		v_ablj=V1212(a,b,l,j);

		vv=v_abil*(2.*v_abjl-v_ablj)+v_abli*(2.*v_ablj-v_abjl);

		epsi_abl = epsi_a + epsi_b - epsi_l;
	
		c_ij += 0.5 * vv
		  *(epsi_abl-epsi_mean)/((epsi_abl-epsi_i)*(epsi_abl-epsi_j));

		f_ij += -0.25 * vv
		  /((epsi_a+epsi_b-epsi_l-epsi_i)*(epsi_a+epsi_b-epsi_l-epsi_j));
	      }
	    }
	  }

	}
      }

      if(c_matrix) c_matrix[m_ind] += c_ij;

      if(f_matrix) {
	f_matrix[XO(i)*symtab->nOcc[inp->sym] + XO(j)] += f_ij;
	if(i!=j)
	  f_matrix[XO(j)*symtab->nOcc[inp->sym] + XO(i)] += f_ij;
      }

      if(inp->debug>9) printf("%4d%4d |%12.8f  %12.8f\n",*i,*j,c_ij,f_ij);
      
      m_ind++;
    }
  }
  if(inp->debug>9) printf("\n#</c11-2>\n\n");
}

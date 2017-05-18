
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_c12_2.c,v 1.2 2003/08/21 09:59:16 joergb Exp $";

#define VV1(K,L,A,J,B,C,S1,S2) \
        ((V1212(K,L,B,C) S1 V1212(K,L,C,B)) \
	 * (V1212(J,A,B,C) S2 V1212(J,A,C,B)))

#define VV2(K,L,A,J,B,I,S1,S2,S3,S4) \
        (S1 V1212(K,I,B,A) * V1212(J,I,B,L) \
	 S2 V1212(K,I,B,A) * V1212(J,I,L,B) \
	 S3 V1212(K,I,A,B) * V1212(J,I,B,L) \
	 S4 2. * V1212(K,I,A,B) * V1212(J,I,L,B))

/* calc_c12_2:
 *
 * calculates the submatrix C12_2 (the 2nd order of the h/2h1p-block) 
 *
 * Spin-free implementation of formula (A11) [JCP 109, 4734 ('98)], see
 * (A3.4),(A3.7) in the PhD thesis of O. Walter. (A3.4b) is negated here.
 *
 * For the ordering of the coefficients in the 2h1p-block see globals.h
 *
 * $Log: calc_c12_2.c,v $
 * Revision 1.2  2003/08/21 09:59:16  joergb
 * added timing statistics in output
 *
 * Revision 1.1  2003/08/19 16:04:19  u85
 * Initial import
 *
 *
 */
void calc_c12_2(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix)
{
  int a_sym,k_sym,l_sym,b_sym,c_sym,i_sym;
  int i_sym_kl,i_sym_lk,kI_sym,kl_sym,lI_sym;
  int *j,*a,*k,*l,*b,*c,*i;
  double epsi_j,epsi_a,epsi_k,epsi_l,epsi_b,epsi_c,epsi_i;
  double *c_ptr1,*c_ptr0=c_matrix;
  double x0,x1,sgn;
  int n,dim_2h1p;

/*   printf("calculating 2nd order of c12 (main/sat-block)\n"); */

  /* this term changes sign in affinity mode */
  sgn = inp->affinity_mode ? -1 : 1 ;

  /*init_time();*/
  dim_2h1p=calc_dim_2h1p(inp,symtab);

  n=0;
  FOR_ALL_2H1P_AKL({
    /* akl part of 2h1p loop (k!=l) */
  
    /*print_time(n,dim_2h1p);*/
  
    c_ptr1 = c_ptr0 + symtab->nOcc[inp->sym];
    
    kl_sym=MULTAB(k_sym,l_sym);
    lI_sym=MULTAB(l_sym,inp->sym);
    
    epsi_a=scf->epsi[*a];
    epsi_k=scf->epsi[*k];
    epsi_l=scf->epsi[*l];
    
    FOR_OCC_IN_SYM_WE(j,inp->sym,epsi_j) {
      x0=x1=0;

      for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	c_sym=MULTAB(b_sym,kl_sym);
	i_sym_kl=MULTAB(b_sym,lI_sym);
	i_sym_lk=MULTAB(b_sym,kI_sym);
	
	FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {
	  
	  /* sum over b and c */
	  FOR_VIR_IN_SYM_WE(c,c_sym,epsi_c) {
	    x0 += SQRT_1_8 * VV1(k,l,a,j,b,c,+,+) / (epsi_k+epsi_l-epsi_b-epsi_c);
	    x1 -= SQRT_3_8 * VV1(k,l,a,j,b,c,-,-) / (epsi_k+epsi_l-epsi_b-epsi_c);
   	  }
	  
	  /* first part of sum over b and i */
	  FOR_OCC_IN_SYM_WE(i,i_sym_kl,epsi_i) {
   	    x0 -= SQRT_1_2 * VV2(k,l,a,j,b,i,+,+,+,-) / (epsi_k+epsi_i-epsi_a-epsi_b);
	    x1 += SQRT_3_2 * VV2(k,l,a,j,b,i,+,-,-,+) / (epsi_k+epsi_i-epsi_a-epsi_b);
	  }
	  
	  /* (l<->k)-part of sum over b and i; Warning i has other sym here */
	  FOR_OCC_IN_SYM_WE(i,i_sym_lk,epsi_i) {
	    x0 -= SQRT_1_2 * VV2(l,k,a,j,b,i,+,+,+,-) / (epsi_l+epsi_i-epsi_a-epsi_b);
	    x1 -= SQRT_3_2 * VV2(l,k,a,j,b,i,+,-,-,+) / (epsi_l+epsi_i-epsi_a-epsi_b);
	  }
	}
      }
      *(c_ptr0++) += sgn * x0;
      *(c_ptr1++) += sgn * x1;
    }
    c_ptr0=c_ptr1;
    n+=2;
  },{
    /* akk part of 2h1p loop */
    
    /*print_time(n,dim_2h1p);*/

    FOR_OCC_IN_SYM_WE(j,inp->sym,epsi_j) {
      x0=0;
      epsi_a=scf->epsi[*a];
      epsi_k=scf->epsi[*k];
      
      for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	c_sym=b_sym;
	i_sym=MULTAB(b_sym,kI_sym);
	
	FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {
	  
	  /* sum over b and c */
	  FOR_VIR_IN_SYM_WE(c,c_sym,epsi_c) {
	    x0 += 0.25 * VV1(k,k,a,j,b,c,+,+) / (epsi_k+epsi_k-epsi_b-epsi_c);
	  }
	  
	  /* sum over b and i */
	  FOR_OCC_IN_SYM_WE(i,i_sym,epsi_i) {
	    x0 -= VV2(k,k,a,j,b,i,+,+,+,-) / (epsi_k+epsi_i-epsi_a-epsi_b);
	  }
	}
      }
      *(c_ptr0++) += sgn * x0;
    }
    n++;
  });

  
/*   if(inp->debug<10) return; */




/*   /\* debug output follows *\/ */

/*   c_ptr0=c_matrix; */
  
/*   fprintf(stdout, */
/* 	  "#<c12-2>\n\n" */
/* 	  "#  a   k   l spin|  j  |  C12(j,akl)\n" */
/* 	  "#----------------+-----+--------------\n"); */
/*   FOR_ALL_2H1P_AKL({ */
/*     FOR_OCC_IN_SYM(j,inp->sym) */
/*       fprintf(stdout,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,1,*j,*(c_ptr0++)); */
/*     FOR_OCC_IN_SYM(j,inp->sym) */
/*       fprintf(stdout,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,2,*j,*(c_ptr0++)); */
/*   },{ */
/*     FOR_OCC_IN_SYM(j,inp->sym) */
/*       fprintf(stdout,"%4d%4d%4d%4d |%4d |%12.8f\n",*a,*k,*l,1,*j,*(c_ptr0++)); */
/*   }); */
/*   fprintf(stdout,"\n#</c12-2>\n\n"); */
}

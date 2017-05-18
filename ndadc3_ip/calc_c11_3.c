
#define SAFE_INT
//#define F_HERM

/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_c11_3.c,v 1.9 2004/01/29 18:48:49 joergb Exp $";

/* calc_c11_3:
 *
 * calculates the submatrix C11_3 (the 3rd order of the hh-block)
 *
 * $Log: calc_c11_3.c,v $
 * Revision 1.9  2004/01/29 18:48:49  joergb
 * corrected formula for hermitian F-matrix (identical to prev. Rel.)
 *
 * Revision 1.7  2003/12/01 14:54:14  joergb
 * F changed into hermitian matrix
 *
 * Revision 1.6  2003/08/21 09:59:16  joergb
 * added timing statistics in output
 *
 * Revision 1.5  2003/08/19 16:50:52  u85
 * small changes in output
 *
 * Revision 1.4  2003/06/12 13:42:03  u85
 * fixed major bug in matrixelements;
 * additional representation of C added
 *
 * Revision 1.3  2003/03/03 14:00:57  u85
 * Just prints matrix in debugmode now
 *
 * Revision 1.2  2003/01/10 16:31:05  u85
 * f_matrix added and bug in c11_3 (h.c. part) fixed
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */
void calc_c11_3(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
		double *c_matrix,double *f_matrix)
{
  int a_sym,b_sym,c_sym,d_sym,l_sym,m_sym,k_sym,ai_sym,bi_sym,ab_sym;
  int *i,*j,*l,*m,*k,*a,*b,*c,*d;
  double epsi_i,epsi_j,epsi_a,epsi_b,epsi_c,epsi_d,epsi_l,epsi_m,epsi_k;
  double epsi_x_epsi,vvv,c_ij,f_ij,f_ji;
  int m_ind,iter;
  
/*   printf("calculating 3rd order of c11 (main/main-block)\n"); */

/*   if(inp->debug>9) { */
/*     fprintf(inp->data_fp,"#<c11-3>\n\n" */
/* 	    "#  i   j | C11(i,j)\n" */
/* 	    "#--------+-----------\n"); */
/*   } */

  m_ind=0;
  iter=symtab->nOcc[inp->sym]*(symtab->nOcc[inp->sym]+1)/2;
  /*init_time();*/

  for(i=symtab->occ[inp->sym];*i!=-1;i++) {
    epsi_i=scf->epsi[*i];

    for(j=symtab->occ[inp->sym];j<=i;j++) {
      /*print_time(m_ind,iter);*/                                                     
      epsi_j=scf->epsi[*j];

      c_ij=0.;
      f_ij=0.;
      f_ji=0.;

      /* C_ij^(A) */
      {
	double v_abil,v_abli,v_cdjl,v_cdlj,v_cdab,v_cdba;

	for(a_sym=1;a_sym<=symtab->nSym;a_sym++) {
	  ai_sym=MULTAB(inp->sym,a_sym);
	  for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	    l_sym=MULTAB(ai_sym,b_sym);
	    ab_sym=MULTAB(a_sym,b_sym);
	    for(c_sym=1;c_sym<=symtab->nSym;c_sym++) {
	      d_sym=MULTAB(ab_sym,c_sym);
	      
	      FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {     
		FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {     
		  FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {     

		    v_abil=V1212(a,b,i,l);
		    v_abli=V1212(a,b,l,i);

		    FOR_VIR_IN_SYM_WE(c,c_sym,epsi_c) {     
		      FOR_VIR_IN_SYM_WE(d,d_sym,epsi_d) {     

			v_cdjl=V1212(c,d,j,l);
			v_cdlj=V1212(c,d,l,j);
			v_cdab=V1212(c,d,a,b);
			v_cdba=V1212(c,d,b,a);
		      
#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 2./ ( 1./( (epsi_a + epsi_b - epsi_l - epsi_i)
				      *(epsi_c + epsi_d - epsi_l - epsi_i))
				 +1./( (epsi_a + epsi_b - epsi_l - epsi_j)
				      *(epsi_c + epsi_d - epsi_l - epsi_j)));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_i)
			  * (epsi_c + epsi_d - epsi_l - epsi_j);
#endif

#ifdef SAFE_INT
			vvv = 0.25 / epsi_x_epsi
  			  * (  2. * v_abil * v_cdjl * v_cdab
			     - 1. * v_abil * v_cdjl * v_cdba
			     - 1. * v_abil * v_cdlj * v_cdab
			     + 2. * v_abil * v_cdlj * v_cdba
			     - 1. * v_abli * v_cdjl * v_cdab
			     + 2. * v_abli * v_cdjl * v_cdba
			     + 2. * v_abli * v_cdlj * v_cdab
			     - 1. * v_abli * v_cdlj * v_cdba);
#else
			vvv = 0.25 / epsi_x_epsi
			  * ( v_abil  * ( v_cdjl  * ( 2.*v_cdab -    v_cdba ) 
					  +v_cdlj * (   -v_cdab + 2.*v_cdba ))
			      +v_abli * ( v_cdjl  * (   -v_cdab + 2.*v_cdba ) 
					  +v_cdlj * ( 2.*v_cdab -    v_cdba )));
#endif

			c_ij += vvv;
			f_ij += vvv / (epsi_c + epsi_d - epsi_l - epsi_i);
			f_ji += vvv / (epsi_a + epsi_b - epsi_l - epsi_j);
		      }
		    }
		  }
		}
	      }
	    
	    }
	  }
	}

      }

      /* C_ij^(B) */
      {
	double v_abil,v_abli,v_acjm,v_acmj,v_lcbm,v_lcmb;

	for(a_sym=1;a_sym<=symtab->nSym;a_sym++) {
	  ai_sym=MULTAB(inp->sym,a_sym);
	  for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	    l_sym=MULTAB(ai_sym,b_sym);
	    for(c_sym=1;c_sym<=symtab->nSym;c_sym++) {
	      m_sym=MULTAB(ai_sym,c_sym);

	      
	      FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {     
		FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {     
		  FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {     

		    v_abil=V1212(a,b,i,l);
		    v_abli=V1212(a,b,l,i);
		  
		    FOR_VIR_IN_SYM_WE(c,c_sym,epsi_c) {     
		      FOR_OCC_IN_SYM_WE(m,m_sym,epsi_m) {     

			v_acjm=V1212(a,c,j,m);
			v_acmj=V1212(a,c,m,j);
			v_lcbm=V1212(l,c,b,m);
			v_lcmb=V1212(l,c,m,b);
		      

#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 2./ ( 1./( (epsi_a + epsi_b - epsi_l - epsi_i)
				      *(epsi_a + epsi_c - epsi_m - epsi_i))
				 +1./( (epsi_a + epsi_b - epsi_l - epsi_j)
				      *(epsi_a + epsi_c - epsi_m - epsi_j)));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_i)
			  * (epsi_a + epsi_c - epsi_m - epsi_j);
#endif

#ifdef SAFE_INT
			vvv = 1. / epsi_x_epsi
			  * (  4. * v_abil * v_acjm * v_lcbm
			     - 2. * v_abil * v_acjm * v_lcmb
			     - 2. * v_abil * v_acmj * v_lcbm
			     + 1. * v_abil * v_acmj * v_lcmb
			     - 2. * v_abli * v_acjm * v_lcbm
			     + 1. * v_abli * v_acjm * v_lcmb
			     + 1. * v_abli * v_acmj * v_lcbm
			     - 2. * v_abli * v_acmj * v_lcmb);
#else
			vvv = 1. / epsi_x_epsi
			  * ( v_abil  * ( v_acjm  * (  4.*v_lcbm - 2.*v_lcmb)
					  +v_acmj * ( -2.*v_lcbm +    v_lcmb))
			      +v_abli * ( v_acjm  * ( -2.*v_lcbm +    v_lcmb)
					  +v_acmj * (     v_lcbm - 2.*v_lcmb)));
#endif

			c_ij += vvv;
			f_ij += vvv / (epsi_a + epsi_c - epsi_m - epsi_i);
			f_ji += vvv / (epsi_a + epsi_b - epsi_l - epsi_j);
		      }
		    }
		  }
		}
	      }

	    }
	  }
	}

      }

      /* C_ij^(C) */
      {
	double v_ablm,v_abml,v_abkj,v_abjk,v_abki,v_abik;
	double v_lmki,v_lmik,v_lmkj,v_lmjk;

	for(a_sym=1;a_sym<=symtab->nSym;a_sym++) {
	  ai_sym=MULTAB(inp->sym,a_sym);
	  for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	    k_sym=MULTAB(ai_sym,b_sym);
	    ab_sym=MULTAB(a_sym,b_sym);
	    for(l_sym=1;l_sym<=symtab->nSym;l_sym++) {
	      m_sym=MULTAB(ab_sym,l_sym);

	      FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {     
		FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {     
		  FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k) {     

		    v_abkj=V1212(a,b,k,j);
		    v_abjk=V1212(a,b,j,k);
		    v_abki=V1212(a,b,k,i);
		    v_abik=V1212(a,b,i,k);

		    FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {     
		      FOR_OCC_IN_SYM_WE(m,m_sym,epsi_m) {     

			v_ablm=V1212(a,b,l,m);
			v_abml=V1212(a,b,m,l);
			v_lmki=V1212(l,m,k,i);
			v_lmik=V1212(l,m,i,k);
			v_lmkj=V1212(l,m,k,j);
			v_lmjk=V1212(l,m,j,k);
		      
#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 1./
			  ( (epsi_a + epsi_b - epsi_k - .5*(epsi_i + epsi_j))
			  / (epsi_a + epsi_b - epsi_l - epsi_m)
			  / (epsi_a + epsi_b - epsi_k - epsi_i)	       
			  / (epsi_a + epsi_b - epsi_k - epsi_j));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_m)
			  * (epsi_a + epsi_b - epsi_k - epsi_j);
#endif

#ifdef SAFE_INT
			vvv = 0.25 / epsi_x_epsi
			  * (  2. * v_ablm * v_abkj * v_lmki
			     - 1. * v_ablm * v_abkj * v_lmik
			     - 1. * v_ablm * v_abjk * v_lmki
			     + 2. * v_ablm * v_abjk * v_lmik
			     - 1. * v_abml * v_abkj * v_lmki
			     + 2. * v_abml * v_abkj * v_lmik
			     + 2. * v_abml * v_abjk * v_lmki
			     - 1. * v_abml * v_abjk * v_lmik);
#else
			vvv = 0.25 / epsi_x_epsi
			  * ( v_ablm  * ( v_abkj  * ( 2.*v_lmki -    v_lmik)
					  +v_abjk * (   -v_lmki + 2.*v_lmik))
			      +v_abml * ( v_abkj  * (   -v_lmki + 2.*v_lmik)
					  +v_abjk * ( 2.*v_lmki -    v_lmik)));
#endif

			c_ij += vvv;
			f_ij += vvv / (epsi_a + epsi_b - epsi_i - epsi_k);
		      

			/* now for the h.c. part */

#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 1./
			  ( (epsi_a + epsi_b - epsi_k - .5*(epsi_i + epsi_j))
			  / (epsi_a + epsi_b - epsi_l - epsi_m)
			  / (epsi_a + epsi_b - epsi_k - epsi_i)
			  / (epsi_a + epsi_b - epsi_k - epsi_j));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_m)
			  * (epsi_a + epsi_b - epsi_k - epsi_i);
#endif

#ifdef SAFE_INT
			vvv = 0.25 / epsi_x_epsi
			  * (  2. * v_ablm * v_abki * v_lmkj
			     - 1. * v_ablm * v_abki * v_lmjk
			     - 1. * v_ablm * v_abik * v_lmkj
			     + 2. * v_ablm * v_abik * v_lmjk
			     - 1. * v_abml * v_abki * v_lmkj
			     + 2. * v_abml * v_abki * v_lmjk
			     + 2. * v_abml * v_abik * v_lmkj
			     - 1. * v_abml * v_abik * v_lmjk);
#else
			vvv = 0.25 / epsi_x_epsi
			  * ( v_ablm  * ( v_abki  * ( 2.*v_lmkj -    v_lmjk)
					  +v_abik * (   -v_lmkj + 2.*v_lmjk))
			      +v_abml * ( v_abki  * (   -v_lmkj + 2.*v_lmjk)
					  +v_abik * ( 2.*v_lmkj -    v_lmjk)));
#endif

			c_ij += vvv;
			f_ji += vvv / (epsi_a + epsi_b - epsi_j - epsi_k);
		      }
		    }
		  }
		}
	      }

	    }
	  }
	}

      }

      /* C_ij^(D) */
      {
	double v_ablm,v_abml,v_bcjm,v_bcmj,v_lcia,v_lcai;
	double v_bcim,v_bcmi,v_lcja,v_lcaj;

	for(a_sym=1;a_sym<=symtab->nSym;a_sym++) {
	  ai_sym=MULTAB(inp->sym,a_sym);
	  for(b_sym=1;b_sym<=symtab->nSym;b_sym++) {
	    bi_sym=MULTAB(inp->sym,b_sym);
	    for(c_sym=1;c_sym<=symtab->nSym;c_sym++) {
	      l_sym=MULTAB(ai_sym,c_sym);
	      m_sym=MULTAB(bi_sym,c_sym);
	      
	      FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b) {     
		FOR_VIR_IN_SYM_WE(c,c_sym,epsi_c) {     
		  FOR_OCC_IN_SYM_WE(m,m_sym,epsi_m) {     

		    v_bcjm=V1212(b,c,j,m);
		    v_bcmj=V1212(b,c,m,j);
		    v_bcim=V1212(b,c,i,m);
		    v_bcmi=V1212(b,c,m,i);
		    
		    FOR_VIR_IN_SYM_WE(a,a_sym,epsi_a) {     
		      FOR_OCC_IN_SYM_WE(l,l_sym,epsi_l) {     

			v_ablm=V1212(a,b,l,m);
			v_abml=V1212(a,b,m,l);
			v_lcia=V1212(l,c,i,a);
			v_lcai=V1212(l,c,a,i);
			v_lcja=V1212(l,c,j,a);
			v_lcaj=V1212(l,c,a,j);
			
#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 1./
			  ( (epsi_b + epsi_c - epsi_m - .5*(epsi_i + epsi_j))
			  / (epsi_a + epsi_b - epsi_l - epsi_m)
			  / (epsi_b + epsi_c - epsi_m - epsi_i)	       
			  / (epsi_b + epsi_c - epsi_m - epsi_j));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_m)
			  * (epsi_b + epsi_c - epsi_m - epsi_j);
#endif

#ifdef SAFE_INT
			vvv = 1. / epsi_x_epsi
			  * (       v_ablm * v_bcjm * v_lcia
			     - 2. * v_ablm * v_bcjm * v_lcai
			     - 2. * v_ablm * v_bcmj * v_lcia
			     + 4. * v_ablm * v_bcmj * v_lcai
			     - 2. * v_abml * v_bcjm * v_lcia
			     + 1. * v_abml * v_bcjm * v_lcai
			     + 1. * v_abml * v_bcmj * v_lcia
			     - 2. * v_abml * v_bcmj * v_lcai);
#else
			vvv = 1. / epsi_x_epsi
			  * ( v_ablm  * ( v_bcjm  * (     v_lcia - 2.*v_lcai )
					  +v_bcmj * ( -2.*v_lcia + 4.*v_lcai ))
			      +v_abml * ( v_bcjm  * ( -2.*v_lcia +    v_lcai )
					  +v_bcmj * (     v_lcia - 2.*v_lcai )));
#endif				

			c_ij += vvv;
			f_ij += vvv / (epsi_b + epsi_c - epsi_i - epsi_m);


			/* now for the h.c. part */
			
#ifdef F_HERM_OLD //FIXME could vanish
			epsi_x_epsi
			  = 1./
			  ( (epsi_b + epsi_c - epsi_m - .5*(epsi_i + epsi_j))
			    / (epsi_a + epsi_b - epsi_l - epsi_m)
			    / (epsi_b + epsi_c - epsi_m - epsi_i)	       
			    / (epsi_b + epsi_c - epsi_m - epsi_j));
#else
			epsi_x_epsi
			  = (epsi_a + epsi_b - epsi_l - epsi_m)
			  * (epsi_b + epsi_c - epsi_m - epsi_i);
#endif

#ifdef SAFE_INT
			vvv = 1. / epsi_x_epsi
			  * (       v_ablm * v_bcim * v_lcja
			     - 2. * v_ablm * v_bcim * v_lcaj
			     - 2. * v_ablm * v_bcmi * v_lcja
			     + 4. * v_ablm * v_bcmi * v_lcaj
			     - 2. * v_abml * v_bcim * v_lcja
			     + 1. * v_abml * v_bcim * v_lcaj
			     + 1. * v_abml * v_bcmi * v_lcja
			     - 2. * v_abml * v_bcmi * v_lcaj);
#else
			vvv = 1. / epsi_x_epsi
			  * ( v_ablm  * ( v_bcim  * (     v_lcja - 2.*v_lcaj )
					  +v_bcmi * ( -2.*v_lcja + 4.*v_lcaj ))
			      +v_abml * ( v_bcim  * ( -2.*v_lcja +    v_lcaj )
					  +v_bcmi * (     v_lcja - 2.*v_lcaj)));
#endif				

			c_ij += vvv;
			f_ji += vvv / (epsi_b + epsi_c - epsi_j - epsi_m);
		      }
		    }
		  }
		}
	      }

	    }
	  }
	}

      }

#ifdef F_HERM
      c_ij += (epsi_i-epsi_j)*(f_ij-f_ji)/2;
#endif

      if(c_matrix) {
	if(!inp->affinity_mode) c_matrix[m_ind] -= c_ij;
	else                    c_matrix[m_ind] += c_ij;
      }

#ifdef F_HERM
      if(f_matrix) {
	f_matrix[XO(i)*symtab->nOcc[inp->sym] + XO(j)] += (f_ij + f_ji) / 2;
	if(i!=j)
	  f_matrix[XO(j)*symtab->nOcc[inp->sym] + XO(i)] += (f_ij + f_ji) / 2;
      }
#else
      if(f_matrix) {
	f_matrix[XO(i)*symtab->nOcc[inp->sym] + XO(j)] += f_ij;
	if(i!=j)
	  f_matrix[XO(j)*symtab->nOcc[inp->sym] + XO(i)] += f_ji;
      }
#endif

      if(inp->debug>9) printf("%4d%4d |%12.8f\n",*i,*j,c_ij);

      m_ind++;
    }
  }
  if(inp->debug>9) printf("\n#</c11-3>\n\n");
}

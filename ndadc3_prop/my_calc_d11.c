/* #include "phis.h" */
#include "my_globals.h"
#include "my_adc_macros.h"
/* #include "calc_indeces.h" */


#define CJ (*c*(*c-1)/2+*j-1)
#define CI (*c*(*c-1)/2+*i-1)

static char *cvsid="$Id: my_calc_d11.c,v 1.0  date time imke Exp$";

/* cap_calc_d11.c:
 *
 * calculates the h/h block of the ISR using ISR code by jb for ISR formulae
 * based on formulae by Trofimov, Schirmer
 * and nd_adc form to build the part of the matrix
 * this part is stored in d11_matrix (must be available for write later on), 'OUTPUT'
 *
 * cap_elements is the matrix containing the (MO/CAP/MO) integrals, INPUT 
 * rho_0 is the the density matrix calculated on beforehand using sigma and q_kl, INPUT 
 * d_null (number!!) is needed for (13a) formula D0(2) ....
 *
 *
 * $Log cap_calc_d11.c, v 1.1 23.04.04 16:52:02 imke
 * 1.1 changed sign before rho_0*d-Term from '+' to '-'
 *     changed indices according to Saschas convention
 *     comment 2 break commands.
 *
 * 1.0 initial version 
 */

/*----------------------------------------------------------------------------------*/
/* NOMENCLATURE NOTE:                                                               */
/*                                                                                  */
/* concerning the formulae: unfortunately they are all written down for i i'        */
/* my nomenclature: keep all indices as they are (even if I may use more int than   */
/* necessary) with two changes:                                                     */
/*                                                                                  */
/* i' -> j                                                                          */
/* such that we need a new name for j -> k                                          */
/*                                                                                  */
/* note: this is a bit different from  jb's nomenclature, which tends to confuse me */
/*----------------------------------------------------------------------------------*/

void cap_calc_d11(struct inp_t *inp, struct scf_t *scf, struct symtab_t *symtab,
                  double *d_null, double *d11_matrix, double *rho_0, double *cap_elements)
{
  int *i, *j, *c, *n, *b, *f, *g, *k;
  int n_sym, b_sym, f_sym, g_sym, k_sym, fi_sym, bi_sym;

  double v_jnfb, v_jnbf, v_fgin, v_fgni;   /* 12a */
  double v_jkfg, v_jkgf;                   /* (with v_fgin, v_fgni)*/
  double v_nkfg, v_nkgf;                   /* (with v_fgin, v_fgni) 12c */
  double v_fgnk, v_fgkn, v_jnfg, v_jngf;   /* (with v_nkfg, v_nkgf ) 12c h.c. */ 

  double epsi_i, epsi_j, epsi_n, epsi_b, epsi_f, epsi_g, epsi_k;

  double sum_1, sum_2;                    /* for intermediate summation, replacing rho_pq, ... */
  double w_ij;                            /* for summation of entry d_ij, cap assumed -> w_ij */

  int dmat_ind=0;                         /* index in d11_matrix */

  int x, index_i, index_j;

  /*------------------------------------------------------------------------------------------- */

  printf("calculating second order one-particle_ISR of CAP for hole/hole block:\n");
  fflush(stdout);

  /*------------------------------------------------------------------------------------------- */


  for(i=symtab->occ[inp->sym];*i!=-1;i++){      /* FOR_OCC_IN_SYM_WE(i,inp->sym,epsi_i){        */
    epsi_i=scf->epsi[*i];
    for(j=symtab->occ[inp->sym];j<=i;j++){     /* FOR_OCC_IN_SYM_WE(j,inp->sym,epsi_j){   j<=i  */
      epsi_j=scf->epsi[*j];
    

      /* initialize */
      w_ij = 0.0;

      /* (13 a) */
      if(*i==*j){
	w_ij += *d_null;
      }

      /* (13b)  */
      w_ij +=  -CAP_PTR(j,i);                 /*d_element(j,i,cap_elements); */


      /* (13c) */
      FOR_VIR_IN_SYM(c,inp->sym){
	
	w_ij -= rho_0[CI] * CAP_PTR(j,c);              /*d_element(j,c,cap_elements);*/
	w_ij -= rho_0[CJ] * CAP_PTR(i,c)           ;   /* h.c. d_element(i,c,cap_elements); */
      }

      /*-- terms composing (13d)  --*/

      /* here all symmetry combinations of b and g have to be included */

      for(b_sym=1; b_sym<=symtab->nSym; b_sym++){
	g_sym=b_sym;
	bi_sym=MULTAB(b_sym,inp->sym);

	/* (12a) */

	FOR_VIR_IN_SYM_WE(b,b_sym,epsi_b){
	  FOR_VIR_IN_SYM_WE(g,g_sym,epsi_g){
	    /*   if(g>b) break;  */                    /* from density ISR where gb matrix is calculated, wrong here */

	    sum_1=0.0;
	    
	    /* loop over all symmteries n=fbi */
	    for(f_sym=1; f_sym<=symtab->nSym; f_sym++){
	      n_sym=MULTAB(f_sym,bi_sym);

	      FOR_VIR_IN_SYM_WE(f,f_sym,epsi_f){
		FOR_OCC_IN_SYM_WE(n,n_sym,epsi_n){

		  v_jnfb=V1212(j,n,f,b);
		  v_jnbf=V1212(j,n,b,f);
		  v_fgin=V1212(f,g,i,n);
		  v_fgni=V1212(f,g,n,i);

		  sum_1 -=
		    ( 2.0*v_jnfb*v_fgin -    v_jnbf*v_fgin
		        - v_jnfb*v_fgni + 2.0*v_jnbf*v_fgni)
		    / (epsi_j + epsi_n - epsi_f - epsi_b)
		    / (epsi_i + epsi_n - epsi_f - epsi_g);
		}
	      }
	    }
	    w_ij += sum_1 * CAP_PTR(b,g);                 /*d_element(b,g,cap_elements);*/
	  }
	}

	/* (12b) */

	/* need: symmtery specifying for n,k */
	k_sym = b_sym;
	n_sym = g_sym;
      
	FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k){
	  FOR_OCC_IN_SYM_WE(n,n_sym,epsi_n){
	    /*if (n>k) break; */                     /* remains from jbs program s.o. */

	    sum_1 = 0.0;

	    /* symmetries: g=fki  = fbi above, USE BI ! */
	    for(f_sym=1; f_sym<=symtab->nSym; f_sym++){
	      g_sym = MULTAB(f_sym,bi_sym);

	      FOR_VIR_IN_SYM_WE(f,f_sym,epsi_f){
		FOR_VIR_IN_SYM_WE(g,g_sym,epsi_g){

		  v_jkfg = V1212(j,k,f,g);
		  v_jkgf = V1212(j,k,g,f);
		  v_fgin = V1212(f,g,i,n);
		  v_fgni = V1212(f,g,n,i);

		  sum_1 += 0.5 *
		    (  2.*v_jkfg*v_fgin -    v_jkgf*v_fgin
		      - v_jkfg*v_fgni + 2.*v_jkgf*v_fgni)
		    / (epsi_j + epsi_k - epsi_f - epsi_g)
		    / (epsi_i + epsi_n - epsi_f -epsi_g);
		}
	      }
	    }
	    w_ij += sum_1 * CAP_PTR(n,k);                    /*d_element(n,k);*/
	  }
	}
      }

      /* (12 c) */
      k_sym = inp->sym;

      FOR_OCC_IN_SYM_WE(k,k_sym,epsi_k){

	sum_1 = 0.0;
	sum_2 = 0.0;

	/* symmtries: n=gfi */
	for(f_sym=1; f_sym<=symtab->nSym; f_sym++){
	  fi_sym = MULTAB(f_sym,inp->sym);
	  for(g_sym=1; g_sym<=symtab->nSym; g_sym++){
	    n_sym = MULTAB(fi_sym,g_sym);

	    FOR_VIR_IN_SYM_WE(f,f_sym,epsi_f){
	      FOR_VIR_IN_SYM_WE(g,g_sym,epsi_g){
		FOR_OCC_IN_SYM_WE(n,n_sym,epsi_n){

		  v_nkfg = V1212(n,k,f,g);
		  v_nkgf = V1212(n,k,g,f);
		  v_fgin = V1212(f,g,i,n);
		  v_fgni = V1212(f,g,n,i);
		  v_fgnk = V1212(f,g,n,k);
		  v_fgkn = V1212(f,g,k,n);
		  v_jnfg = V1212(j,n,f,g);
		  v_jngf = V1212(j,n,g,f);

		  /* !h.c. */
		  sum_1 -= 0.25*
		    ( + v_nkfg*v_fgin - 2.0*v_nkgf*v_fgin
		      - 2.0*v_nkfg*v_fgni + v_nkgf*v_fgni)
		    / (epsi_n + epsi_k - epsi_f - epsi_g)
		    / (epsi_i + epsi_n - epsi_f - epsi_g);

		  sum_2 -= 0.25*
		    ( + v_fgnk*v_jnfg - 2.0*v_fgkn*v_jnfg
		      -2.0*v_fgnk*v_jngf + v_fgkn*v_jngf)
		    / (epsi_n + epsi_k - epsi_f -epsi_g)
		    / (epsi_j + epsi_n - epsi_f -epsi_g);
		}
	      }
	    }
	  }
	}
	w_ij += sum_1 * CAP_PTR(j,k) + sum_2 * CAP_PTR(i,k); /* d_element(j,k,cap_elements)   / d_element(i,k,cap_elements) */
      }
    

    if(d11_matrix) d11_matrix[dmat_ind] += w_ij; 

      dmat_ind++;
    }
  }


/*   /\* debug output *\/ */
/*   if(inp->debug > 1){ */

/*     for(x = 0; x < dmat_ind; x++){ */
/*       calc_indexpair_sp(x,&index_i,&index_j); */
/*       fprintf(inp->data_fp, "global index of h/h-matrix: %i ij_index: (%i,%i), Element: %16.14f\n", x, index_i, index_j, d11_matrix[x]); */
/*     } */
/*   } */


}






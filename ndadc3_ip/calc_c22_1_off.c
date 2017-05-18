
/* #include "phis.h" */
#include "globals.h"
#include <math.h>

static char *cvsid="$Id: calc_c22_1_off.c,v 1.2 2003/09/02 11:11:23 joergb Exp $";

/* this implements \delta_KM * V * V needed severeal times */
#define DELTA_V_TERM(K,M,A,N,B,L,S00,S01,S10,S11) \
        if(K==M) { \
	  v1 = V1212(A,N,B,L); \
	  v2 = V1212(A,N,L,B); \
	  x00 += S00 prefact*(-v1 + 0.5 * v2); \
	  x01 -= S01 prefact*(SQRT_3_4 * v2); \
	  x10 -= S10 prefact*(SQRT_3_4 * v2); \
	  x11 += S11 prefact*(-v1 + 1.5 * v2); \
	}



/* calc_c22_1_off: (see below)
 *
 * calculates the submatrix C22_1 (the 1st order of the 2h1p/2h1p-block) 
 *
 * Spin-free implementation of formula (A12) [JCP 109, 4734 ('98)], see
 * (A3.5),(A3.8) in the PhD thesis of O. Walter. (A3.5c) is negated here.
 *
 * For the ordering of the coefficients in the 2h1p-block see globals.h
 *
 * $Log: calc_c22_1_off.c,v $
 * Revision 1.2  2003/09/02 11:11:23  joergb
 * debug output in calc_c22_1_cols corrected
 *
 * Revision 1.1  2003/08/19 16:04:19  u85
 * Initial import
 *
 *
 */



/* calc_c22_1_cols:
 *
 * calculates one (if m==n) or two columns (spin 1 and 2) of C22_1 up to the 
 * actual row number given in 'row' (since the matrix is symmetric, only the 
 * upper part is calculated)
 */
void calc_c22_1_cols(struct inp_t *inp,struct scf_t *scf,
		     struct symtab_t *symtab,
		     int *m,int *n,int *b,int row,double *c_col0,double *c_col1)
{
  int k_sym,l_sym,a_sym,kI_sym;
  int *k,*l,*a;
  double v_mnkl,v_mnlk,v1,v2;
  double *c_ptr0=c_col0,*c_ptr1=c_col1;
  double prefact,sgn;

  double x00,x01,x10,x11;

  /* this term changes sign in affinity mode */
  sgn = inp->affinity_mode ? -1 : 1 ;

  /* if m==n we multiply all terms by sqrt(1/2), see O. Walter (A3.8) */
  prefact=1.;
  if(*m==*n) prefact=SQRT_1_2;

  /* now loop through the column 
   * (warning: the loop breaks with a GOTO, when the diagonal is reached
   *           => we just calulate the upper part of the matrix here)
   */
  FOR_ALL_2H1P_AKL({
    /* akl part of 2h1p loop (k!=l) */

    /* FIXME: should be done outside of loop */    
    c_ptr0[0] = 0; c_ptr0[1] = 0; c_ptr1[0] = 0; c_ptr1[1] = 0;

    x00=x01=x10=x11=0;
	    
    if(a==b) {
      v_mnkl=V1212(m,n,k,l);
      v_mnlk=V1212(m,n,l,k);
    
      x00 += prefact*(v_mnkl + v_mnlk);
      x11 += prefact*(v_mnkl - v_mnlk);
    }
	    
    DELTA_V_TERM(k,m,a,n,b,l,+,+,+,+);
    DELTA_V_TERM(l,n,a,m,b,k,+,-,-,+);
    DELTA_V_TERM(l,m,a,n,b,k,+,-,+,-);
    DELTA_V_TERM(k,n,a,m,b,l,+,+,-,-);

    *(c_ptr0++) += sgn * x00;
    *(c_ptr0++) += sgn * x01;
    *(c_ptr1++) += sgn * x10;
    *(c_ptr1++) += sgn * x11;

    if(c_ptr0 - c_col0 > row) goto debug_and_return;  /* oh-oh goto! */
  },{
    /* akk part of 2h1p loop */

    c_ptr0[0] = 0; c_ptr1[0] = 0; /* FIXME: should done outside of loop */

    x00=x10=x01=x11=0;
	    
    if(a==b) 
      x00 += prefact*2.*V1212(m,n,k,k);

    DELTA_V_TERM(k,m,a,n,b,k,+,+,+,+);
    DELTA_V_TERM(k,n,a,m,b,k,+,-,-,+);
    DELTA_V_TERM(k,m,a,n,b,k,+,-,+,-);
    DELTA_V_TERM(k,n,a,m,b,k,+,+,-,-);

    x00 *= SQRT_1_2;
    x10 *= SQRT_1_2;
    
    *(c_ptr0++) += sgn * x00;
    *(c_ptr1++) += sgn * x10;

    if(c_ptr0 - c_col0 > row) goto debug_and_return;  /* oh-oh goto! */
  });




 debug_and_return:
  
  if(inp->debug<10) return;




/*   /\* debug output follows *\/ */

/*   c_ptr0=c_col0; */
/*   c_ptr1=c_col1; */

/*   /\* FIXME: printing complete (not packed) spin matrix on diagonal *\/ */
/*   if(*m!=*n) { */
/*     FOR_ALL_2H1P_AKL({ */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,1,*(c_ptr0++)); */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,2,*(c_ptr0++)); */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,2,1,*(c_ptr1++)); */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,2,2,*(c_ptr1++)); */

/*       if(c_ptr0 - c_col0 > row) return; */
/*     },{ */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,1,*(c_ptr0++)); */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,2,1,*(c_ptr1++)); */

/*       if(c_ptr0 - c_col0 > row) return; */
/*     }); */
/*   } */
/*   else { */
/*     FOR_ALL_2H1P_AKL({ */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,1,*(c_ptr0++)); */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,2,*(c_ptr0++)); */

/*       if(c_ptr0 - c_col0 > row) return; */
/*     },{ */
/*       fprintf(inp->data_fp,"%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/* 	      *b,*m,*n,*a,*k,*l,1,1,*(c_ptr0++)); */

/*       if(c_ptr0 - c_col0 > row-1) return; */
/*     }); */
/*   } */
}



/* calc_c22_1_off :
 *
 * calculates the submatrix C22_1 (the 1st order of the 2h1p/2h1p-block) 
 * using the column subroutine above
 */
void calc_c22_1_off(struct inp_t *inp,struct scf_t *scf,
		    struct symtab_t *symtab,int dim_2h1p,double *c_22)
{
  int k_sym,l_sym,a_sym,kI_sym;
  int *k,*l,*a;
  double *c_col0,*c_col1,*c_ptr=c_22;
  int row,i;

  printf("calculating 1st order off-diagonal elements of c22 (sat/sat-block)\n");

  /* allocate space for two columns of C22 */
  c_col0=(double*)malloc((dim_2h1p+1)*sizeof(double));
  c_col1=(double*)malloc((dim_2h1p+1)*sizeof(double));


/*   /\* start debug output *\/ */
/*   if(inp->debug>=10)  */
/*     fprintf(inp->data_fp, */
/* 	    "#<c22-1-off>\n\n" */
/* 	    "# a'  k'  l' |   a   k   l |Spins|  C22(a'k'l',akl)\n" */
/* 	    "#------------+-------------+-----+-------------------\n"); */


  /* first column has one (!) (diagonal-) element */
  row=1;

  /* loop over rows and calculate the columns */
  FOR_ALL_2H1P_AKL({
    /* akl part of 2h1p loop (k!=l) */

    calc_c22_1_cols(inp,scf,symtab,k,l,a,row,c_col0,c_col1);
    
    for(i=0;i<row;i++) *(c_ptr++) = c_col0[i]; row++;
    for(i=0;i<row;i++) *(c_ptr++) = c_col1[i]; row++; 
  },{
    /* akk part of 2h1p loop */

    calc_c22_1_cols(inp,scf,symtab,k,l,a,row,c_col0,c_col1);
  
    for(i=0;i<row;i++) *(c_ptr++) = c_col0[i]; row++;
  });


  /* end of debug output */
/*   if(inp->debug>=10)  */
/*     fprintf(inp->data_fp,"\n#</c22-1-off>\n\n"); */

  free(c_col0);free(c_col1);
}



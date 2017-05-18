/* #include "phis.h" */
#include "my_globals.h"
#include <math.h>
/* #include "my_adc_macros.h" */

static char *cvsid="$Id: cap_calc_d22_off.c,v 1.1 22.04.04 18:26:00 imke Exp$";

/* cap_calc_c22_off: (not yet implemented)
 *
 * calculates the off diagonal part of the ISR of a one-Particle operator
 * 
 * 
 * spin-free formulae by imke
 * 
 * ordering of coefficients and "how to go through matrix" similar to nd_adc(3) by joergb! 
 *
 * $Log: cap_calc_d22_off.c,v $
 * 
 */

/* cap_calc_d22_off_cols:
 * 
 * calculates two colums (spin 1 and 2 ) of off-diagonal part of 2h1p/2h1p-block up to the
 * actual row number given in 'row' 
 *
 * note: only one of the two cols is really used if the col. treated belongs to k==l case 
 * (here named n==m) !!
 *
 * m,n,b: akl of col to be calculated 
 * cap_elements (MO/CAP/MO)
 *
 * d_col1, d_col0: columns to be calculated 
 *
 * $Log: cap_calc_d22_off_cols.c,v $
 * 1.1 22.04.04 18:26:00 imke
 *     exchange contributions due to c and e for consistency with Saschas formulae
 *
 * 1.0 initial version
 */


void cap_calc_d22_off_cols(struct inp_t *inp, struct scf_t *scf, struct symtab_t *symtab,
                           int *m, int *n, int *b, int row, double *d_col0, double *d_col1, 
                           double *cap_elements)
{
  int k_sym, l_sym, a_sym, kI_sym;
  int *k, *l, *a;

  double *d_ptr0=d_col0, *d_ptr1=d_col1;

  /*double sgn;*/ /* not yet used, for affinity mode ... to be implemented ... */

  double x00, x01, x10, x11;

  /*---------------------------------------------------------------------*/
  /* loop through the column bmn:                                        */
  /* this loop uses jb's goto ...                                        */
  /* to calculate only upper part of the matrix ...                      */
  /*---------------------------------------------------------------------*/

  FOR_ALL_2H1P_AKL({

    /* akl part of 2h1p loop within column (k!=l)                  */
    /* case m!=n:                                                  */
    /* |1 0|            |-1 0|                                     */
    /* |0 1| (b,c,e)    | 0 1|  ((d,f)                             */
    /* -> two cols used !!                                         */ 
    /*                                                             */
    /* case m==n:                                                  */
    /* | 1/sqrt(2)|          |- 1/sqrt(2)|                         */
    /* | 1/sqrt(6)| (b,c,e)  |  1/sqrt(6)|  (d,f)                  */
    /* -> two cols calculated but only one used by calling routine */
                            


    d_ptr0[0] = 0; d_ptr0[1] = 0; d_ptr1[0]=0; d_ptr1[1]=0;
    x00=x01=x10=x11=0.0;

    /* here: m!=n and m==n are only different py prefactor -> calc everything and check for pref. then */
    /* TODO: replace the if's by else if's because if one is fulfilled non of the other should be poss. */
    /* or better by switch ? */

    /* (1a) only in diagonal, nothing to be done here */

    /* (1b) */
    if(*k==*m && *l==*n){
      x00 += CAP_PTR(a,b);             /*d_element(a,b,cap_elements);*/
      x11 += CAP_PTR(a,b);             /*d_element(a,b,cap_elements);*/
    }

    /* (1c,d,e,f) */
    if(*a==*b){

      /* c */
      if(*k==*m){
	x00 -= CAP_PTR(l,n);           /*d_element(l,n,cap_elements);*/
	x11 -= CAP_PTR(l,n);           /*d_element(l,n,cap_elements);*/
      }
       
      /* e */
      if(*l==*n){
	x00 -= CAP_PTR(k,m);           /*d_element(k,m,cap_elements);*/
	x11 -= CAP_PTR(k,m);           /*d_element(k,m,cap_elements);*/
      }

      /* d */
      if(*k==*n){
	x00 += -CAP_PTR(m,l);          /*d_element(m,l,cap_elements);*/
	x11 +=  CAP_PTR(m,l);          /*d_element(m,l,cap_elements);*/
      }

      /* f */
      if(*l==*m){
	x00 += -CAP_PTR(k,n);          /*d_element(k,n,cap_elements);*/
	x11 +=  CAP_PTR(k,n);          /*d_element(k,n,cap_elements);*/
      }
    }
    
    /* m == n ? */
    if(*m==*n){
      x00 = x00 * SQRT_1_2;
      x01 = x11 * SQRT_1_6; /* attention: the "minus" and "plus" SHOULD be correct like this, but ... */
         /* note: this does only work if col0 is used by calling routine, else chaos!*/ 
    }

    *(d_ptr0++) += x00;
    *(d_ptr0++) += x01;
    *(d_ptr1++) += x10;
    *(d_ptr1++) += x11;

    if(d_ptr0 - d_col0 > row) goto debug_and_return;    /* goto from jb, if possible replace by st nicer*/

                  
  },{
    /* k==l part of a,k,l loop                                               */
    /* -> only one entry to each of the probably two colums!                 */
    /*                                                                       */
    /* case m!=n :                                                           */
    /* |1/sqrt(2)  1/sqrt(6)|  (b,c,e)  |-1/sqrt(2) 1/sqrt(6)| (d,f)         */
    /*                                                                       */
    /* case m==n:                                                            */
    /*  1    (b,c,e)           0 (d,f)                                       */

    /* here: first if about m==n             ....                            */

    d_ptr0[0] = 0; d_ptr1[0] = 0;
    x00=x10=0.0;

    if(*m!=*n){

      /* a, nothing to by done */

      /* TODO: b, there should never be an contribution ...? */

      if(*k==*m && *l==*n){
	x00 += SQRT_1_2 * CAP_PTR(a,b);         /*d_element(a,b,cap_elements);*/
	x10 += SQRT_1_6 * CAP_PTR(a,b);         /*d_element(a,b,cap_elements);*/
      }

      /* c,e,d,f */
      if(*a==*b){


	/* c */
	if(*k==*m){
	  x00 -= SQRT_1_2 * CAP_PTR(l,n);       /*d_element(l,n,cap_elements);*/
	  x10 -= SQRT_1_6 * CAP_PTR(l,n);       /*d_element(l,n,cap_elements);*/
	}
	
	/* e */
	if(*l==*n){
	  x00 -= SQRT_1_2 * CAP_PTR(k,m);       /*d_element(k,m,cap_elements);*/
	  x10 -= SQRT_1_6 * CAP_PTR(k,m);       /*d_element(k,m,cap_elements);*/
	}

	/* d */
	if(*k==*n){
	  x00 += -SQRT_1_2 * CAP_PTR(l,m);      /*d_element(l,m,cap_elements);*/
	  x10 +=  SQRT_1_6 * CAP_PTR(l,m);      /*d_element(l,m,cap_elements);*/
	}

	/* f */
	if(*l==*m){
	  x00 += -SQRT_1_2 * CAP_PTR(k,n);      /*d_element(k,n,cap_elements);*/
	  x10 +=  SQRT_1_6 * CAP_PTR(k,n);      /*d_element(k,n,cap_elements);*/
	}
      }
    }
    else
      {

	/* a, nothing */
 
	/* b */
	if(*k==*m && *l==*n){
	  x00 += CAP_PTR(a,b);                 /*d_element(a,b,cap_elements);*/
	  x10 += CAP_PTR(a,b);                 /*d_element(a,b,cap_elements);*/
	}

	/* c, d, e, f */

	if(a==b){

	  /* c */
	  if(*k==*m){
	    x00 -= CAP_PTR(l,n);               /*d_element(l,n,cap_elements);*/
	    x10 -= CAP_PTR(l,n);               /*d_element(l,n,cap_elements);*/
	  }

	  /* e */
	  if(*l==*n){
	    x00 -= CAP_PTR(k,m);               /*d_element(k,m,cap_elements);*/
	    x10 -= CAP_PTR(k,m);               /*d_element(k,m,cap_elements);*/
	  }

	  /* d */
	  /*if(*k==*n){
	  x00 += 0; 
	  x10 += 0; 
	  }*/

	  /* f */
	  /*if(*l==*m){
	  x00 += 0;
	  x10 += 0;
	  }*/
	}
      }

    *(d_ptr0++) += x00;
    *(d_ptr1++) += x10;  /* second col is calculated here, but should not be */
                         /* in use because this is the k==l part of the */
		         /* calling routine */

    if(d_ptr0 - d_col0 > row) goto debug_and_return; /* goto */
  });

 debug_and_return:
  ;
/*   if(inp->debug<10) return; */

/*   /\* debug *\/ */

/*   d_ptr0=d_col0; */
/*   d_ptr1=d_col1; */

/*   printf("debugging off-disg of 2h1p/2h1p ...... \n"); */
  
/*   if(*m!=*n) { */
/*     FOR_ALL_2H1P_AKL({ */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,1,*(d_ptr0++)); */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,2,*(d_ptr0++)); */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,2,1,*(d_ptr1++)); */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,2,2,*(d_ptr1++)); */

/*       if(d_ptr0 - d_col0 > row) return; */
/*     },{ */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,1,*(d_ptr0++)); */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,2,1,*(d_ptr1++)); */
      
/*       if(d_ptr0 - d_col0 > row) return; */
/*     }); */
/*   } */
/*   else { */
/*     FOR_ALL_2H1P_AKL({ */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,1,*(d_ptr0++)); */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,2,*(d_ptr0++)); */

/*       if(d_ptr0 - d_col0 > row) return; */
/*     },{ */
/*       printf("%4d%4d%4d |%4d%4d%4d |%2d%2d | %12.8f\n", */
/*               *b,*m,*n,*a,*k,*l,1,1,*(d_ptr0++)); */

/*       if(d_ptr0 - d_col0 > row-1) return; */
/*     }); */
/*   } */



}



#define  PI   3.14159265358979312
#define  BEGIN {
#define  END   };
#include <stdio.h>
#include <math.h>
/*#include "fortran.h"*/
/*---------------------------------------------------------------------------*/
/*                                                         UR   15.01.95     */
/*---------------------------------------------------------------------------*/

void intgauss_ (double *wert, double *alpha, int *n )

/*---------------------------------------------------------------------------*/
/*
/*     Routine for the calculation of the integral of a Gaussian :
/*              
/*     wert = _/^  X^n  e^{ -alpha * X^2 } 
/*
/*---------------------------------------------------------------------------*/

{

     int    k, i;
     double a;
     double x[12] = { 1.0, 0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875 };

     if ( *n % 2 == 1 ) 
        *wert = 0;
     else
   {  
        k = *n / 2;
        a = 1;
        for ( i=0; i<k; ++i) a *= *alpha;
        *wert = x[k] * sqrt(PI / *alpha) / a ;
   }

     

}










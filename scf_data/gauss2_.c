#define  BEGIN {
#define  END   };
#include <stdio.h>
#include <math.h>
/*#include "fortran.h"*/
/*---------------------------------------------------------------------------*/
/*                                                         UR   01.12.94     */
/*---------------------------------------------------------------------------*/

  void gauss2_ (double *ao, double *alpha, 
               double *bo, double *beta, 
               double *co, double *expo, double *faktor ) 

/*---------------------------------------------------------------------------*/
/*
/*     Routine for the calculation of a new Gaussian :
/*
/*     faktor * Gauss(co,exo) = Gauss(ao,alpha) * Gauss(bo,beta)
/*
/*---------------------------------------------------------------------------*/

{

     *expo   = *alpha + *beta;

     *co     = ( *ao * *alpha + *bo * *beta ) / *expo;

     *faktor = *alpha * *beta  * (*bo - *ao) * (*bo - *ao) / *expo;

     *faktor = exp( -(*faktor) );

     

}










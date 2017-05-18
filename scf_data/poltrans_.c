#define  BEGIN {
#define  END   };
#include <stdio.h>
#include <math.h>
/*#include "fortran.h"*/
/*---------------------------------------------------------------------------*/
/*                                                         UR   15.01.95     */
/*---------------------------------------------------------------------------*/

  void poltrans_ (double f[], double *alpha, double *beta, int *n )

/*---------------------------------------------------------------------------*/
/*
/*     Routine for the transformation of an polynomial:
/*              
/*     {sum} a_i ( x - a )^i   ->   {sum} b_i ( x - b )^i
/* 
/*---------------------------------------------------------------------------*/

{

     int    i, j, k;
     double wert, diff;
  
     diff = *beta - *alpha;
     k    = 1;

     for ( j=0; j<*n+1; j++)
     BEGIN
     wert = 0;
     for ( i=*n; i>j-1; i--) 
     BEGIN
        wert = wert * diff + f[i];
        f[i] = (i-j) * f[i]; 
     END
     f[j] = wert  / k;
     k    = (j+1) * k;
     END;

     

}










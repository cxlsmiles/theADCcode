#define  BEGIN {
#define  END   };
#include <stdio.h>
#include <math.h>
/*#include "fortran.h"*/
/*---------------------------------------------------------------------------*/
/*                                                         UR   15.01.95     */
/*---------------------------------------------------------------------------*/

void polmulti_ (double f[], double g[], int *n, 
                             double h[], int *m )
/*---------------------------------------------------------------------------*/
/*
/*     Routine for the transformation of an polynomial:
/*              
/*     {sum} a_i x^i  *  {sum} b_i x^i  ->  {sum} c_i x^i
/* 
/*---------------------------------------------------------------------------*/

{

     int    i, j;

     for ( i=0; i<*n+*m+1; f[i++] = 0 );
     for ( i=0; i<*n+1; i++ )
     for ( j=0; j<*m+1; j++ )
     f[i+j] = f[i+j] + g[i] * h[j];

     

}










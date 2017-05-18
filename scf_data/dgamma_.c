#include <stdio.h>
#include <math.h>
double  dgamma_(double *xt)
{
  double dgamma, lg;
  extern int signgam;
   
  lg=gamma(*xt);
    dgamma=signgam*exp(lg); 

	     return dgamma;
	     }

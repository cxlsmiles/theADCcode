#include <stdlib.h>

#include <stdio.h>
#include "jbutils.h"

static char *cvsid="$Id: make_matrix.c,v 1.3 2003/08/19 16:54:01 u85 Exp $";

/* make_matrix:
 *
 * generates matrix of given dimension with some additional features
 *
 * parameters:
 *   int rows:        number of rows
 *   int cols:        number of columns
 *   char mode:       matrix mode ('g'=general, 'p'=symmetric packed)
 *   char fill:       fill matrix ('n'=no, 'z'=zeroes)
 *   double **matrix: returns matrix
 *
 * $Log: make_matrix.c,v $
 * Revision 1.3  2003/08/19 16:54:01  u85
 * can create non quadratic matrices (e.g. c12) now
 *
 * Revision 1.2  2003/01/10 16:31:06  u85
 * f_matrix added and bug in c11_3 (h.c. part) fixed
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */
void make_matrix(int rows,int cols,char mode,char fill,double **matrix)
{
  int len,i;

  switch(mode) {
  case 'g': len=rows*cols;     break;
  case 'p': 
    if(rows!=cols) { 
      fflush(stdout);
      fprintf(stderr,"make_matrix: matrix must be squared for packed mode\n"); 
      exit(42);  
    }
    len=rows*(rows+1)/2; 
    break;
  default: 
    fflush(stdout);
    fprintf(stderr,"make_matrix: mode %c not implemented\n",mode);
    exit(42);
  }

  /* allocate matrix */
  *matrix = (double*)malloc_we(len*sizeof(double));

  /* fill with zeroes */
  if(fill=='z') {
    for(i=0;i<len;i++)
      (*matrix)[i]=0.;
  }
  else if(fill!='n') {
    fflush(stdout);
    fprintf(stderr,"make_matrix: fill-mode %c not implemented\n",mode);
    exit(42);
  }
}

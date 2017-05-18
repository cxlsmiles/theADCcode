#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include "calc_indeces.h"

static char *cvsid="$Id: my_calc_indeces.c,v 1.0 2004/04/27 16:18:14 imke Exp $";


/* my_calc_indeces.c:
 *
 * functions for calculating local index-pairs (i,j) beginning with i=j=1
 * from global c-matrix-index beginning at 0
 *
 * matrix storage supported:
 * 
 * a) calc_indexpair_sp(int glob_index, int *i, int *j)
 *   upper triangular by column of a squared matrix
 *
 * b) calc_indexpair_cf(int num_rows, int glob_index, int *i, int *j)
 *   full matrix by column , number_of_rows has to be input 
 *
 *
 * $Log: my_calc_indeces.c,v $
 * Revision 1.0 2004/04/27 16:18:14 imke
 * initial version
 */


void calc_indexpair_sp(int glob_index, int *i, int *j)
{

  /* calculates the indices (i,j) of a matrix element of a upper triagonal
   * by row matrix beginning with (1,1) from global c-index beginning with 0
   *
   * estimates j from glob_index making use of 
   * glob_index <= j*(j+1)/2 -1
   * -> trial_index_float=-0.5+sqrt(0.25+2*(glob_index+1))
   *    trial_index_int=int(trail_index_float) => j-1 if non-diag j if diag
   *    set trial_index_int=trial_index_int+1
   *
   * calculate i from glob_index=j*(j-1)/2-1 +i
   * -> trial_index_i=glob_index + 1 -(trial_index_int*(trail_index_int-1))/2
   * 
   * case diagonal: j==0: set i=j=trial_index_int-1
   */ 

  double float_help;
  int index_help1;
  int index_help2;

  float_help=-0.5+sqrt(0.25+2.0*(glob_index+1));
  index_help1=(int)float_help +1;
  index_help2=glob_index +1 -(index_help1*(index_help1-1))/2;
  if(index_help2!=0){
    *j=index_help1;
    *i=index_help2;
  }
  else{
    *i=*j=index_help1-1;
  }
}

void calc_indexpair_cf(int num_rows, int glob_index, int *i, int *j)
{
  /* calculates indexpair (i,j) of a full matrix with num_rows rows
   * stored by col. from global index 
   *
   * column_index-1 = int((float)global_index:(float)(num_rows))
   * row_index -1 = (int)global_index % (int) num_rows
   */

  *j= (int)((float)glob_index/num_rows) +1;
  *i= (glob_index % num_rows) +1;
}



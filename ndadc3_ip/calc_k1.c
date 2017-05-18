
#include "globals.h"

static char *cvsid="$Id: calc_k1.c,v 1.5 2003/08/19 16:50:52 u85 Exp $";

/* calc_k1:
 *
 * calculates the submatrix K1 (the 0th order of the hh-block) and the 
 * 0th order of the hh-block of f1
 *
 * $Log: calc_k1.c,v $
 * Revision 1.5  2003/08/19 16:50:52  u85
 * small changes in output
 *
 * Revision 1.4  2003/03/03 14:22:42  u85
 * error corrected in c_matrix
 *
 * Revision 1.3  2003/01/31 15:52:11  u85
 * debug output removed
 *
 * Revision 1.2  2003/01/10 16:31:03  u85
 * f_matrix added and bug in c11_3 (h.c. part) fixed
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */
void calc_k1(struct inp_t *inp,struct scf_t *scf,struct symtab_t *symtab,
	     double *c_matrix,double *f_matrix)
{
  int i;

/*   printf("calculating 0th order of c11 (main/main-block)\n"); */

  for(i=0;i<symtab->nOcc[inp->sym];i++) {
    if(c_matrix) c_matrix[i*(i+1)/2+i] -= scf->epsi[symtab->occ[inp->sym][i]];
    if(f_matrix) f_matrix[i*symtab->nOcc[inp->sym]+i] += 1;
  }
}

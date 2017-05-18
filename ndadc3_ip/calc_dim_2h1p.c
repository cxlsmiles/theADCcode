
/* #include "phis.h" */
#include "globals.h"

static char *cvsid="$Id: calc_dim_2h1p.c,v 1.1.1.1 2002/12/04 15:17:17 u85 Exp $";

/* calc_dim_2h1p:
 *
 * calculates the number of 2h1p-states
 *
 * $Log: calc_dim_2h1p.c,v $
 * Revision 1.1.1.1  2002/12/04 15:17:17  u85
 * Imported sources
 *
 */
int calc_dim_2h1p(struct inp_t *inp,struct symtab_t *symtab)
{
  int sym_i,sym_j,dim_2h1p;

  /* calculate the number of 2h1p elements */
  dim_2h1p=0;
  for(sym_i=1;sym_i<=symtab->nSym;sym_i++) 
    for(sym_j=1;sym_j<sym_i;sym_j++) 
      dim_2h1p += 
	2 * symtab->nOcc[sym_i]
	* symtab->nOcc[sym_j]
	* symtab->nVir[MULTAB(MULTAB(sym_i,sym_j),inp->sym)];

  for(sym_i=1;sym_i<=symtab->nSym;sym_i++)
    dim_2h1p +=
      2 * symtab->nOcc[sym_i] * (symtab->nOcc[sym_i]-1) / 2
      * symtab->nVir[inp->sym]
      + symtab->nOcc[sym_i] * symtab->nVir[inp->sym];

  return(dim_2h1p);
}
 

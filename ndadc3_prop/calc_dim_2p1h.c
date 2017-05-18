
#include "phis.h"
#include "my_globals.h"

static char *cvsid="$Id: calc_dim_2p1h.c,v 1.1.1.1 2002/12/04 15:17:17 u85 Exp $";

/* calc_dim_2p1h:
 *
 * calculates the number of 2p1h-states
 *
 * $Log: calc_dim_2p1h.c,v $
 * Revision 1.1.1.1  2002/12/04 15:17:17  u85
 * Imported sources
 *
 */
int calc_dim_2p1h(struct inp_t *inp,struct symtab_t *symtab)
{
  int sym_i,sym_j,dim_2p1h;

  /* calculate the number of 2p1h elements */
  dim_2p1h=0;
  for(sym_i=1;sym_i<=symtab->nSym;sym_i++) 
    for(sym_j=1;sym_j<sym_i;sym_j++) 
      dim_2p1h += 
	2 * symtab->nVir[sym_i]
	* symtab->nVir[sym_j]
	* symtab->nOcc[MULTAB(MULTAB(sym_i,sym_j),inp->sym)];
  
  for(sym_i=1;sym_i<=symtab->nSym;sym_i++)
    dim_2p1h +=
      2 * symtab->nVir[sym_i] * (symtab->nVir[sym_i]-1) / 2
      * symtab->nOcc[inp->sym]
      + symtab->nVir[sym_i] * symtab->nOcc[inp->sym];

  return(dim_2p1h);
}
 

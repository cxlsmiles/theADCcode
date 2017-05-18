/* $Id: make_symtab.h,v 1.2 2003/12/01 14:54:55 joergb Exp $
 *
 * $Log: make_symtab.h,v $
 * Revision 1.2  2003/12/01 14:54:55  joergb
 * parameters changed (make_symtab does not need the
 * whole inp structure => usable in other programs)
 *
 * Revision 1.1  2003/10/05 08:49:46  joergb
 * splitted into globals.h, adc_macros.h, get_input.h,
 *
 */



/* threshold for orbital occupancy */
#define SCF_THRESH 1e-6


#include "read_scf_data.h"

/* symtab_t:
 * structure for symmetry information
 */
struct symtab_t {
  int    nSym;            /* number of symmetries */
  int    nOcc[9];         /* number of occupied orbitals per sym */
  int    nVir[9];         /* number of virtual orbitals per sym */
  int    *occ[9];         /* symmetry mapping for occupied orbitals per sym */
  int    *vir[9];         /* symmetry mapping for virtual orbitals per sym */
};


/* void make_symtab(int sym,int affinity_mode,int debug,FILE *data_fp, */
/* 		 struct scf_t *scf,struct symtab_t *symtab); */



#include "./guk.h"
/*#include "typedefs.h"*/

/* sort_gukscf.c:
 * 
 * The scf-vectors of GUK are sorted in a strange way. Here they are 
 * reordered. 
 *
 * $Log: sort_gukscf.c,v $
 * Revision 1.4  2003/04/24 14:05:49  joergb
 *  * guk_aux.c (guk_getpara): New function, determines GUK parameters
 *    from actual dfile
 *    (guk_sizeof_subsection): New function, likewise
 *  * configure.ac: Likewise (getpara.pl call removed)
 *  * guk.h: Likewise
 *  * guk_init.c: Likewise
 *  * phis_guk_ao.c: Likewise
 *  * phis_guk_geometry.c: Likewise
 *  * typeinits.c: New file, initializes the GUK typdefs using the
 *    determined GUK parameters
 *  * typedefs.h: Likewise
 *  * Makefile.in: Likewise
 *
 *  * guk_init.c: now detects dfile generated with GUK5 (needed for
 *    information in section type 2
 *  * guk5_read_type2.c: New file, reads several matrices in section
 *    type 2 from GUK5 dfiles: S,T,H,x,y,z
 *  * phis_guk_overlap.c: Likewise
 *
 * Revision 1.3  2002/11/15 20:02:22  alext
 * * guk.h (phis_guk_geometry): New function.
 * (nAtoms): New global variable.
 * * guk_init.c: Determines nAtoms.
 * * phis_guk_info.c: Likewise.
 * * phis_guk_geometry.c: New file.
 * * Makefile.in: Likewise.
 * * phis_guk_overlap.c: Error return changed: -nBas*nBas -> -nBas
 *
 * * phis_guk_ao.c: fixed another bug considering
 * f-functions in lists "l_tab" and "f_tab".
 * * sort_gukscf.c: bug in sorting with n>1 fixed.
 *
 * Revision 1.2  2002/04/10 22:25:58  alext
 * Latest patches from Joerg
 *
 * Revision 1.1  2002/03/24 22:20:03  alext
 * * guk_aux.c (guk_get_section): Hardcoded length of section type 3,
 * needs to be fixed.
 * * phis_guk_info.c: Now returns the correct number of atoms.
 *
 * * sort_gukscf.c: New file.
 * * Makefile.in: New files added.
 * * phis_get_scfvec.c: Call sort_gukscf() to fix AO ordering in SCF
 * vectors.
 *
 * Revision 1.2  2002/02/06 15:03:28  u85
 * comments added/corrected
 *
 * Revision 1.1  2002/02/01 14:04:46  u85
 * C files added
 *
 */

void sort_gukscf(double *gukscf,Type_3 *gukvec)
{
  int i,j,k,l;
  double *vector;


  /* allocating temporary vector */

  vector=(double*)malloc(gukvec->n_gtos*sizeof(double));
  if(vector==NULL){
    fprintf(stderr,"malloc: Error allocating %d B\n",gukvec->n_gtos*sizeof(double));
    perror("malloc:");
    exit(1);
  }


  for(i=0;i<gukvec->n_gtos;i++) {
    for(j=0;j<gukvec->n_gtos;j++) 
      vector[j]=0;


    /* transformation for symmetry equivalent orbitals */
    for(k=0,j=0;j<gukvec->n_gtos;j++) {
      for(l=0;l<gukvec->ntran[j];l++,k++) {
	vector[gukvec->itran[k]-1]+=
	  gukvec->ctran[k]*gukscf[i*gukvec->n_gtos + j];
      }
    }

    /* The sorted vector is copied to the pos. of the orig. vector */
    for(j=0;j<gukvec->n_gtos;j++) 
      gukscf[i*gukvec->n_gtos+j]=vector[j];
  }

  free(vector);
}

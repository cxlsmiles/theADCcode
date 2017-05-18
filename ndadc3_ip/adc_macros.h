/* $Id: adc_macros.h,v 1.2 2004/01/29 17:04:39 joergb Exp $
 *
 * $Log: adc_macros.h,v $
 * Revision 1.2  2004/01/29 17:04:39  joergb
 * corrected XO and XV to work in all symetries
 *
 * Revision 1.1  2003/10/05 08:49:46  joergb
 * splitted into globals.h, adc_macros.h, get_input.h,
 *
 */


#include <stdio.h>


/* constants */
#define SQRT_1_2 0.707106781186547462
#define SQRT_1_8 0.353553390593273731
#define SQRT_3_2 1.22474487139158894
#define SQRT_3_4 0.866025403784438597
#define SQRT_3_8 0.61237243569579447
#define SQRT_2   1.41421356237309515


/* Use these macros to "trick" Joerg's code 
   that phis has been loaded in the normal way*/
#define MULTAB(A,B) Multab_((A),(B))
/* we like the physical 1212 mapping of integrals */
#define V1212(A,B,C,D) vpqrs_((A),(C),(B),(D))

/* loops over all occupied/virtual orbitals in given symmetry B
 * (orbital number inside the loop is *A)
 */
#define FOR_OCC_IN_SYM(A,B) for((A)=symtab->occ[(B)];*(A)!=-1;(A)++)
#define FOR_VIR_IN_SYM(A,B) for((A)=symtab->vir[(B)];*(A)!=-1;(A)++)


/* same as above, with energy, the energy of *A is returned in C */
#define FOR_OCC_IN_SYM_WE(A,B,C) \
     for((A)=symtab->occ[(B)];*(A)!=-1 && (((C)=scf->epsi[*(A)]) || 1);(A)++)
#define FOR_VIR_IN_SYM_WE(A,B,C) \
     for((A)=symtab->vir[(B)];*(A)!=-1 && (((C)=scf->epsi[*(A)]) || 1);(A)++)


/* FOR_ALL_2H1P_AKL(A,B):
 *
 * loop over the 2h1p-indices of the matrix and execute A; or B; for each
 * index-tripel a,k,l:  A, if k!=l - B, if k==l
 *
 * 1.) loop over states with l_sym<k_sym (=> k!=l) => A;
 * 2.) loop over l_sym==k_sym, l<k => A;
 * 3.) loop over k=l (=> k_sym==l_sym) => B;
 *
 * uses the following variables:
 *   int k_sym    symmetry of occ. orb. k
 *   int l_sym    symmetry of occ. orb. l
 *   int a_sym    symmetry of vir. orb. a
 *   int kI_sym   symmetry k * sym of ADC-matrix 
 *
 *   int *k       pointer to occ. orb. k
 *   int *l       pointer to occ. orb. l
 *   int *a       pointer to vir. orb. a
 *
 * (the macro should be readable and tells you the full truth)
 */
#define FOR_ALL_2H1P_AKL(A,B) { \
  for(k_sym=1;k_sym<=symtab->nSym;k_sym++) { \
    kI_sym=MULTAB(inp->sym,k_sym); \
    for(l_sym=1;l_sym<k_sym;l_sym++) { \
      a_sym=MULTAB(kI_sym,l_sym); \
\
      FOR_OCC_IN_SYM(k,k_sym) { \
	FOR_OCC_IN_SYM(l,l_sym) { \
	  FOR_VIR_IN_SYM(a,a_sym) { A; }}}} \
    l_sym=k_sym; \
    a_sym=inp->sym; \
\
    FOR_OCC_IN_SYM(k,k_sym) { \
      for(l=symtab->occ[l_sym];l<k;l++) { \
	FOR_VIR_IN_SYM(a,a_sym) { A; }} \
\
      FOR_VIR_IN_SYM(a,a_sym) { B; }}} \
}


/* X?(A) = orbital number of A in this symmetry */
#define XO(A) ((A)-symtab->occ[scf->sym[*(A)]])
#define XV(A) ((A)-symtab->vir[scf->sym[*(A)]]+symtab->nOcc[scf->sym[*(A)]])




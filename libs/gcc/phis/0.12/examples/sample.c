/** AT ****************************************************************
 *                                                                    *
 * This little sample application calculates the MP2 correction to    *
 * the ground state energy of a closed-shell system.                  *
 *                                                                    *
 * It serves three purposes:                                          *
 * - to demonstrate how to use the PHIS library,                      *
 * - to quickly check new backends of this library, and               *
 * - to give beginners in the field of quantum chemistry a feeling of *
 *   how a typical calculation works.                                 *
 *                                                                    *
 * The same thing in FORTRAN can be found in fsample.f.               *
 *                                                                    *
 **************************************************************** AT **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "phis.h"

#define MAXBAS     1024
#define MAXSYM     9
#define OPTIONS    "mg"

void print_usage(void);

int main(int argc, char **argv)
{

/**********************************************************************/
/* declarations */
	long cap;
	int opt, backend, len;
	
	int nSym, nBas, nAtoms;
	double E_hf;         /* HF energy */
	double *epsi;        /* orbital energies, commonly denoted by epsilon */
	double *occ;         /* occupation numbers */
	int *sym;            /* symmetry labels */

	int *i, *j;          /* indices for hole orbitals */
	int *a, *b;          /* indices for particle orbitals */
	int p;               /* index for arbitrary orbitals */

	int iSym, jSym;      /* the corresponding symmetry labels */
	int aSym, bSym;
	int pSym;

	double iajb, ibja;   /* two-particle integrals */
	double MP2;          /* final MP2 energy */

	int ixj;             /* direct product of symmetry labels */

	/* I don't want to clobber this demo with technical stuff too much,
	   so I use statically allocated arrays here. */
	int list_of_occ[MAXSYM][MAXBAS], list_of_vir[MAXSYM][MAXBAS];

/**********************************************************************/
/* initializations */

	if (argc < 2) print_usage();
	while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
		switch (opt) {
		case 'm' :
			backend = MOLCAS |(SYM_BLOCKED << 16);
			break;
		case 'g' :
			backend = 1;
			break;
		default:
			print_usage();
		}
	}
	cap = phis_init(&backend);
	/* Check if all information is available */
/*	if (~(cap | ~(HAVE_INFO | HAVE_EPSI | HAVE_SYM | HAVE_OCC))) {
		fprintf(stderr, "Sorry, the backend is incomplete. Terminating.\n");
		exit(-1);
	}
*/
	phis_get_info(&nSym, &nBas, &nAtoms);
	printf("no. of irreps: %i\n", nSym);
	printf("no. of basis functions: %i\n", nBas);

	/* IMPORTANT
	 * We need arrays that are one element too large, since C arrays
	 * count from 0, but orbitals are labelled starting with 1. */
	len = nBas; 
	if ((epsi = (double *) malloc((len+1)*sizeof(double))) == NULL ||
	    (sym = (int *) malloc((len+1)*sizeof(double))) == NULL ||
	    (occ = (double *) malloc((len+1)*sizeof(double))) == NULL) {
		perror("malloc");
		exit(-1);
	}

	/* Here we pass the adress of the second element, 
	   which makes e.g. epsi[1] refer to the first orbital. */
	phis_get_epsi(&E_hf, epsi+1, &len);
	phis_get_sym(sym+1, &len);
	phis_get_occ(occ+1, &len);

	if (len < 0 ) {
		/* This should never happen, but you never know... */
		printf("Oops, something is wrong with the number of basis functions.\n");
		printf("The value changed from %i to %i. I will continue with the latter value.\n",
		       nBas, -len);
		nBas = -len;
	}

	printf("\nno.\tsym\tocc\torbital energy\n");
	printf("-----------------------------------------------\n");
	for (p = 1; p <= nBas; p++) {
		printf("%i\t%i\t%2.1f\t%f\n", p, sym[p], occ[p], epsi[p]);
	}
	printf("\nelectronic HF-SCF energy: %f\n", E_hf);

/*
 * You will need a list of occupied (hole) and virtual (particle) orbitals
 * to loop over. You could avoid this if you order your orbitals in a 
 * special way so that orbitals of one kind are aligned in memory, but
 * this makes the code harder to read and much less flexible.
 */
	{
		/* these are local arrays needed nowhere else */ 
		int nOcc[MAXSYM] = { 0 }, nVir[MAXSYM] = { 0 };

		for (p = 1; p <= nBas; p++) {
			pSym = sym[p];
			if (occ[p] == 2.0) {
				list_of_occ[pSym][nOcc[pSym]] = p;	
				nOcc[pSym]++;
			} else if (occ[p] == 0.0) {
				list_of_vir[pSym][nVir[pSym]] = p;
				nVir[pSym]++;
			} else {
				fprintf(stderr, "Sorry, this is only for closed-shell systems.\n");
				fprintf(stderr, "orbital %i had occupation %f\n", p, occ[p]);
				exit(-1);
			}
		}
		
		/* mark the end of the lists */
		for (pSym = 1; pSym <= nSym; pSym++) {
			list_of_occ[pSym][nOcc[pSym]] = -1;
			list_of_vir[pSym][nVir[pSym]] = -1;
		}
	}

/*
 * If you are new to C you might wonder how the loop over orbitals below works:
 * Note that i is a pointer, so initially i points to list_of_occ[pSym][0],
 * the loop terminates if i points to the value -1 (if i is a pointer,
 * *i is the value pointed to), and finally i++ makes i point to the 
 * next element in list_of_occ[pSym][]
 */
	printf("\norbitals by occupation and symmetry:\n");
	for (pSym = 1; pSym <= nSym; pSym++) {
		printf("irrep = %i\n", pSym);
		printf("occupied: ");
		for (i = list_of_occ[pSym]; *i != -1; i++)
			printf(" %i ", *i);
		printf("\nvirtual: ");
		for (a = list_of_vir[pSym]; *a != -1; a++)
			printf(" %i ", *a);
		printf("\n\n");
	}

/**********************************************************************/
	
	MP2 = 0;
	phis_load_Vpqrs();  /* Reads 2e integrals into memory, so we can use Vpqrs() */

	/* outer loop: irreps
	 * the whole thing must be totally symmetric, so bSym is fixed
	 */
	for (iSym = 1; iSym <= nSym; iSym++) {
		for (jSym = 1; jSym <= nSym; jSym++) {
			ixj = MULTAB(iSym,jSym);
			for (aSym = 1; aSym <= nSym; aSym++) {
				bSym = MULTAB(ixj,aSym); 

				/* inner loop over orbitals */
				for (i = list_of_occ[iSym]; *i != -1; i++) {
					for (j = list_of_occ[jSym]; *j != -1; j++) {

						for (a = list_of_vir[aSym]; *a != -1; a++) {
							for (b = list_of_vir[bSym]; *b != -1; b++) {

								iajb = Vpqrs(i,a,j,b);
								ibja = Vpqrs(i,b,j,a);

								MP2 += (2*iajb*iajb - iajb*ibja)
									/ (epsi[*i]+epsi[*j]-epsi[*a]-epsi[*b]);
							}
						}
					}
				}
			}
		}
	}
	printf("MP2 correction: %20.12f\n", MP2 );
	exit(0);
}

void
print_usage(void) 
{
	printf("\n\
sample [-m|-g]\n\
-m\t selects MOLCAS backend\n\
-g\t selects GAMESS-UK backend\n" );
	
	exit(-1);
}

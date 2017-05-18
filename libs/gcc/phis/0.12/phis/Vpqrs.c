#include "config.h"
#include "phis_guk.h"

#ifndef DEBUG
# define NDEBUG        /* Ignore assertions unless in debug mode. */
#else
static int *SymDbg;    /* Remember symmetry labels for debugging. */
#endif

#include <assert.h>
int        nInt;
static int *iv;
static int *iu;
static int *iw;
static double *TwoE;

/*
 * local functions
 */

/* Initialize the tables for integral lookup */
int ixset(int *nBas, int *sym);

/* Return the canonical index for i,j,k,l */
inline int canonicalize(const int *i, const int *j, const int *k, const int *l);


int 
ixset(int *nBasTot, int *SymLab)
/*- AT -----------------------------------------------------------------*
 *                                                                      *
 * Sets up the index fields iu(), iv() and iw() for integral lookup,    *
 * See, for example, I. Shavitt in "Methods of Electronic Structure     *
 * Theory", Henry F. Schaefer III (Ed.)                                 *
 *                                                                      *
 *----------------------------------------------------------------- AT -*/
{

	int i, ii, j, jx, k, lx, l, lm;

	/* this trick is for FORTRAN compatibility */
	int n = *nBasTot;
	int *sym = SymLab - 1;	/* now sym[1] points to SymLab[0] */

	int ix = 0;
//	int nInt = 0;

	int ntt2;
        nInt=0;

	if (iv != NULL) {
		free(iv);
		free(iu);
		free(iw);
	}

	ntt2 = n * (n + 1) / 2;
	if ((iv = (int *) malloc((n + 1) * sizeof(int))) == NULL ||
	    (iu = (int *) malloc((ntt2 + 1) * sizeof(int))) == NULL ||
	    (iw = (int *) malloc((ntt2 + 1) * sizeof(int))) == NULL) {
		fprintf(stderr, "Allocation error.\n");
		exit(-1);
	}

	for (i = 1; i <= n; i++) {
		iv[i] = ix;
		ii = sym[i];
		for (j = 1; j <= i; j++) {
			jx = MULTAB(sym[j], ii);
			ix = ix + 1;
			iu[ix] = nInt;
			for (k = 1; k <= i; k++) {
				lx = MULTAB(sym[k], jx);
				lm = k;
				if (k == i)
					lm = j;
				for (l = 1; l <= lm; l++) {
					if (sym[l] == lx)
						nInt = nInt + 1;
				}
			}
			iw[ix] = nInt - iu[ix];
		}
	}
#ifdef _DEBUG_
        printf("no. of integrals:       %i\n", nInt);
#endif         
	return (nInt);
}


inline int 
canonicalize(const int *i, const int *j, const int *k, const int *l)
{
	int p, q, r, s, pqrs;
	int x, y;

	if (*i >= *j) {
		p = *i;
		q = *j;
	} else {
		p = *j;
		q = *i;
	}

	if (*k >= *l) {
		r = *k;
		s = *l;
	} else {
		r = *l;
		s = *k;
	}

	x = iv[p] + q;
	y = iv[r] + s;
	if (x >= y) { 
		pqrs = iu[x] + iw[y];
	} else {
		pqrs = iu[y] + iw[x];
	}

	return (pqrs);
}


void
guk_phis_load_Vpqrs(void)
{
	int nSym, nOrb, nAtoms;
	int *SymLab = NULL;

	guk_phis_get_info(&nSym, &nOrb, &nAtoms);

	/* get the symmetry labels */
	if ((SymLab = (int *) malloc(nOrb * sizeof(int))) == NULL)  exit(-2);
	guk_phis_get_sym(SymLab, &nOrb);

	nInt = ixset(&nOrb, SymLab);

	if ((TwoE = (double *) malloc((nInt+1) * sizeof(double))) == NULL)
		exit(-2);

	/* Fill the array of 2e integrals */
/* VVP : cycle without checking && conditions */
	for (;;) {
		int p, q, r, s, pqrs;
		double x;
		guk_phis_get_next_Vpqrs(&p, &q, &r, &s, &x);
		if (s == 0) break;
		pqrs = canonicalize(&p, &q, &r, &s);
		TwoE[pqrs] = x; 		
	}
#ifdef DEBUG
	SymDbg = SymLab - 1;
#else
	free(SymLab);
#endif
	return;
}


double 
guk_Vpqrs(const int *i, const int *j, const int *k, const int *l)
{
	int pqrs;

#ifdef DEBUG
	{
		int Vpqrs_sym =
			MULTAB(
				MULTAB(SymDbg[*i],SymDbg[*j]),
				MULTAB(SymDbg[*k],SymDbg[*l]));

		if (Vpqrs_sym != 1) {
			fprintf(stderr, 
				"Vpqrs: non-symmetric integral requested: "
				"(%i %i | %i %i)\n", *i, *j, *k, *l);
			exit(-1);
		}
	}
#endif

	pqrs = canonicalize(i,j,k,l);

	return TwoE[pqrs];
}

double
guk_Vordered(const int *i, const int *j, const int *k, const int *l)
{
	int x, y, pqrs;

	assert(*i >= *j);
	assert(*k >= *k);

	x = iv[*i] + *j;
	y = iv[*k] + *l;

	assert(x >= y);
	pqrs = iu[x] + iw[y];

	return TwoE[pqrs];
}


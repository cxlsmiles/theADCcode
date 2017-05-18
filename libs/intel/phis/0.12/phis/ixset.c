#include <phis_guk.h>

static int *iv;
static int *iu;
static int *iw;

int ixset(int *nBasTot, int *SymLab)
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
	int nInt = 0;

	int ntt2;

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

	return (nInt);
}

int canonicalize(const int *i, const int *j, const int *k, const int *l)
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

	if (p < r || (p == r && q < s)) {
		int t;
		t = p; p = r; r = t;
		t = q; q = s; s = t;
	}

	x = iv[p] + q;
	y = iv[r] + s;
	pqrs = iu[x] + iw[y];

	return (pqrs);
}


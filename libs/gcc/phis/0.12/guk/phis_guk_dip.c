#include <string.h>

#include "phis.h"
#include "./guk.h"

static void transform(double *D, double *T, int nMO, int nAO);

void
phis_guk_dip(double *X, double *Y, double *Z, int *n_max)
{
	int nGtos = section3->n_gtos;
	int bufsize = (nGtos+1)*nGtos/2;

	double *buf, *T;

	/* make sure the matrices are large enough. */
	if (*n_max < nAct) {
		*n_max = -nAct;
		return;
	}

	buf = (double *) malloc_we(bufsize*sizeof(double));
	T = (double *) malloc_we(nGtos*nGtos*sizeof(double));

	guk_read_type2(buf, bufsize, DIPX);
	copy_and_expand(T, buf, nGtos);
	transform(X, T, nAct, nGtos);

	guk_read_type2(buf, bufsize, DIPY);
	copy_and_expand(T, buf, nGtos);
	transform(Y, T, nAct, nGtos);

	guk_read_type2(buf, bufsize, DIPZ);
	copy_and_expand(T, buf, nGtos);
	transform(Z, T, nAct, nGtos);

	free(buf); free(T);
	*n_max = nAct;
	return;
}


static void
transform(double *M, double *A, int nMO, int nAO) 
{
/* Transforms the matrix A from the AO basis to the MO basis
 * and returns the result in M.
 * FIXME: this assumes that there is no stride in M!
 */
	char opA = 'N', opB = 'N';
	double one = 1.0, zero = 0.0;
	double *C, *T;

	/* get the SCF vectors */
	C = (double *) malloc_we(nAO*nMO*sizeof(double));
	phis_guk_scfvec(C, &nMO, &nAO);

	T = (double *) malloc_we(nAO*nAO*sizeof(double));

	/* T = A.C */
	dgemm(&opA, &opB, &nAO, &nMO, &nAO, &one,
	      A, &nAO, C, &nAO, &zero, T, &nAO);

	/* M = Transpose(C).T */
	opA = 'T';
	dgemm(&opA, &opB, &nMO, &nMO, &nAO, &one,
	      C, &nAO, T, &nAO, &zero, M, &nMO);

	free(T); free(C);
	return;
}

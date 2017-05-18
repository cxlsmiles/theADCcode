#include <string.h>
#include "./guk.h"
#include "phis.h"

void
phis_guk_overlap(double *S, int *n_max)
{
	int nGtos = section3->n_gtos;
	int bufsize = (nGtos+1)*nGtos/2;

	double *buf;

	/* make sure we have enough space in S */
	if (*n_max < nGtos) {
		*n_max = -nGtos;
		return;
	}

	buf = (double *) malloc_we (bufsize*sizeof(double));
	guk_read_type2(buf, bufsize, OLAP);
	copy_and_expand(S, buf, nGtos);

	free(buf);
	*n_max = nGtos;
	return;
}

#include "./guk.h"

void
phis_guk_occ(double *occ, int *n)
{
	int i;

	if (*n < nAct) {
		*n = -nAct;
		return;
	}

	for (i=0; i<nAct; i++) 
		occ[i] = section3->occ_num[map_active[i]];
	*n = nAct;

	return;
}

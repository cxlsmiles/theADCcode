#include "./guk.h"

void
phis_guk_epsi(double *E_hf, double *epsi, int *n)
{
	int i;
	double E_nuc;

	if (*n < nAct) {
		*n = -nAct;
		return;
	}

	guk_read_type2(&E_nuc, 1, ENUC);

	for (i=0; i<nAct; i++) 
		epsi[i] = section3->orb_eng[map_active[i]];
	*n = nAct;
	*E_hf = section3->ehf_tot - E_nuc;

	return;
}

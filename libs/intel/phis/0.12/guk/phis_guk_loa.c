#include "./guk.h"

void
phis_guk_loa(int *list, int *n)
{
	int i;
	
	if (*n < nAct) {
		*n = -nAct;
		return;
	}

	for (i=0; i<nAct; i++)                 /* in the application,     */
		list[i] = map_active[i] + 1;   /* orbitals count from one */

	*n = nAct;
	return;
}

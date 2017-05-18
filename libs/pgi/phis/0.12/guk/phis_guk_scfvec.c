#include <string.h>
#include "./guk.h"

void
phis_guk_scfvec(double *C, int *n_max, int *stride)
{
	int i, nGtos=section3->n_gtos;

	Type_3 *not_needed;

	int n2 = nBas * nGtos;

	FILE *dump_fp;
	double *ptr, *buf;

	/* some preparations */
	if ((dump_fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		perror( "Cannot open dumpfile" ); 
		exit( 1 );
	}
	
	/* first check boundaries */
	if (*n_max < nAct) 
		*n_max = -nAct;

	if (*stride < nGtos)
		*stride = -nGtos;

	if ( *stride<0 || *n_max<0 )
		return;


	/* better positioning by first reading the vectors section */
	not_needed = guk_get_section(dump_fp,3);
	free(not_needed->data);
	free(not_needed);


    	if ((buf = (double *) malloc(n2*sizeof(double))) == NULL)
		exit(0);

	if (guk_read_data(dump_fp, (char *) buf, 0, n2) != n2) {
		printf("Something wrong reading SCF vectors.\n");
		exit(1); 
	}

	fclose( dump_fp );

	/* fix up individual SCF vectors */
	sort_gukscf(buf,section3);

	/* take care of frozen and deleted orbitals */
	ptr = C;
	for (i=0; i<nAct; i++) {
		memcpy(ptr, buf+map_active[i]*nGtos, nGtos*sizeof(double));
		ptr += nGtos;
	}
	free(buf);

	*n_max = nAct;
	*stride = nGtos;
	return;
}

#include "./guk.h"
#include <math.h>

void
phis_guk_sym(int *sym, int *n)
{
	int i, j;

	int *map_sym;
	int map_sym_ident[]={0,1,2,3,4,5,6,7,8};
	int map_sym4[]={0,1,3,4,2};
	int map_sym8[]={0,1,8,7,2,6,3,4,5};

	double *scfvectors;
	Type_3 *not_needed;
	int n2;
	FILE *dump_fp;

	if (*n < nAct) {
		*n = -nAct;
		return;
	}

	/* reading scf vectors (needed for determining symmetry from AOs) */
	if ((dump_fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		perror( "Cannot open dumpfile" ); 
		exit( 1 );
	}
  
	not_needed = guk_get_section(dump_fp,3);
	free(not_needed->data);
	free(not_needed);

	n2=nBas*nBas;
	if ((scfvectors = (double *) malloc(n2*sizeof(double))) == NULL)
		exit(0);
  
	if (guk_read_data(dump_fp, (char *) scfvectors, 0, n2) != n2) {
		printf("Something wrong reading SCF vectors.\n");
		exit(1); 
	}

	fclose( dump_fp );


	/* do symmetry mapping from GUK to PHIS symmetries */
	switch(section51->n_irrep) {
	case 8:  map_sym=map_sym8; break;
	case 4:  map_sym=map_sym4; break;
	default: map_sym=map_sym_ident;
	}

  
	/* read the MO symmetries */
	for (i=0; i<nAct; i++) 
		sym[i] = map_sym[section51->symlabels_mo[map_active[i]]];


	/* get symmetry using the SCF vectors and the AO symmetries */
	for (i=0; i<nAct; i++) {
		if(sym[i] == 0) {
			for (j=0; j<nBas; j++) {
				if(fabs(scfvectors[map_active[i]*nBas+j]) > 1.0e-4) {
					if(!sym[i]) 
					  {sym[i]=map_sym[section51->symlabels_ao[j]];
					    
					  }
					else if (sym[i]!=map_sym[section51->symlabels_ao[j]]) {
						fprintf(stderr,
							"Inconsistency in the AO symmetry labels detected!\n"
							"(If you use 'HARMONIC ON', please do not use 'BYPASS HF' in the integral transformation)\n"
							"Aborting.\n");
						exit(-2);
					}
				}
			}
			if(!sym[i]) {
				fprintf(stderr,
					"Unable to determine the symmetry of orbital %d\n",i+1);
				exit(-2);
			}
		}
	}


	*n = nAct;

	free(scfvectors);

	return;
}

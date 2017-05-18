#include <math.h>

#include "phis.h"
#include "./guk.h"

int nBas, nAct, nAtoms;
int *map_active;

Type_3 *section3;
Type_51 *section51;

int guk_mxprim;
int guk_maxaqm;
int guk_maxat;
int guk_mxshel;
int guk_maxorb;

static int query_guk_version(FILE *fp);

int
guk_init_hook(int flags) 
{
	int i;
	FILE *dump_fp;

	guk_interface.info = &phis_guk_info;
	guk_interface.epsi = &phis_guk_epsi;
	guk_interface.sym = &phis_guk_sym;
	guk_interface.occ = &phis_guk_occ;
	guk_interface.next_Vpqrs = &phis_guk_next_Vpqrs;
	guk_interface.scfvec = &phis_guk_scfvec;
	guk_interface.overlap = &phis_guk_overlap;
	guk_interface.geometry = &phis_guk_geometry;
	guk_interface.ao = &phis_guk_ao;
	guk_interface.list_active = &phis_guk_loa;
	guk_interface.dip = &phis_guk_dip;

	if (flags & SYM_BLOCKED) {
		guk_interface.cap = (HAVE_INFO | HAVE_EPSI | HAVE_SYM | HAVE_OCC
				 | HAVE_SCFVEC | HAVE_LOA
				 | HAVE_GEOMETRY | HAVE_DIP);
	} else {
		guk_interface.cap |=
			(HAVE_INFO | HAVE_EPSI | HAVE_SYM | HAVE_OCC | HAVE_NEXTINT 
			 | HAVE_SCFVEC | HAVE_OVERLAP | HAVE_AO | HAVE_LOA 
			 | HAVE_GEOMETRY | HAVE_DIP);
	}

	if ((dump_fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		fprintf(stderr, "Cannot open dumpfile %s\n", dumpfile_name); 
		exit(1);
	}

/*
 * initialize some data
 */

	/* get parameters maxorb,mxshel,... from dfile */
	guk_getpara(dump_fp);

	if (section51) {free(section51->data);free(section51);}
	if ((section51 = (Type_51 *)guk_get_section(dump_fp, SYMMETRY_DATA)) == NULL) {
		perror("guk_get_section");
		exit(1); 
	}

	if (section3) {free(section3->data);free(section3);}
	if ((section3 = (Type_3 *)guk_get_section(dump_fp, SCF_DATA)) == NULL) {
		perror("guk_get_section");
		exit(1); 
	}
	nBas = section3->n_vec;

	/*** Set up the mapping from external (active) orbitals to the SCF basis. ***/
	{

		Type_1005 *section1005;
		
		if ((section1005 = (Type_1005 *)guk_get_section(dump_fp, 1005)) == NULL) {
			fprintf(stderr, "You need to specify RUNTYPE TRANSFORM in your GAMESS run.\n");
			perror("guk_get_Section");
			exit(1); 
		}

		nAct = section1005->nAct;

		if (map_active) free(map_active);
		map_active = (int *) malloc_we(nAct*sizeof(int));

		for (i=0; i<nAct; i++) 
			map_active[i] = section1005->list_of_act[i] - 1;  /* arrays count from zero internally */

		/* If symmetry blocking is requested we need an additional reordering. */
		if (flags & SYM_BLOCKED) {
			int *tmp, *map;
			int n = nAct;
			int j,k;

			tmp = (int *) malloc_we(nAct*sizeof(int));
			map = (int *) malloc_we(nAct*sizeof(int));
			phis_guk_sym(tmp, &n);
			if (n < 0) {
				fprintf(stderr, "guk_init: problem reading symmetries.\n");
				exit(1);
			}

			/* Generate sym->active mapping. */
			for (i=1, k=0; i<=section51->n_irrep; i++)
				for (j=0; j<nAct; j++)
					if (tmp[j] == i)
						map[k++] = j;

			/* Concatenate mappings sym->act and act->bas. */
			for (i=0; i<nAct; i++) 
				tmp[i] = map_active[map[i]];

			/* Save the result. */
			for (i=0; i<nAct; i++)
				map_active[i] = tmp[i];

			free(map); free(tmp);
		}
		
		free(section1005->data);
		free(section1005);
	}

	/* GUK version dependent settings */
	if (query_guk_version(dump_fp) == 5) {
		printf("Warning: dfile seems to generated using GUK version 5\n"
		       " => entering GUK-5 compability mode\n");
		guk_read_type2 = &guk5_read_type2;
	} else
		guk_read_type2 = &guk6_read_type2;

	
	/* read number of atoms from section type 1 */
	{
		Type_1 *section1;

		section1 = guk_get_section(dump_fp,1);
		nAtoms=section1->n_atoms;
		free(section1->data);
		free(section1);
	}

	fclose(dump_fp);

	return guk_interface.cap;
}

int
query_guk_version(FILE *fp)
{
	double *ptr = (double *) Guk_buffer;
	double thresh = 1e-6;

	int sec = guk_sec_by_type(CORE_AO);
	int pos = Summary.entry[sec].pos;

	(void) guk_get_block(fp, pos);

	/* if beginning of overlap matrix follows, this seems to be Guk version 5 */
	if(fabs(ptr[4]-1.)<thresh && fabs(ptr[6]-1.)<thresh && fabs(ptr[9]-1.)<thresh)
		return 5;
	else
		return 6;
}

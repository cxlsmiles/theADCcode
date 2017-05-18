#include "config.h"
#include "phis_guk.h"

bool check( char* );

/*
 * some global things
 */

/* group multiplication table, symmetry labels count from one */
int F77_FUNC(multabguk,MULTABGUK)[64] =
{
	1, 2, 3, 4, 5, 6, 7, 8,
	2, 1, 4, 3, 6, 5, 8, 7,
	3, 4, 1, 2, 7, 8, 5, 6,
	4, 3, 2, 1, 8, 7, 6, 5,
	5, 6, 7, 8, 1, 2, 3, 4,
	6, 5, 8, 7, 2, 1, 4, 3,
	7, 8, 5, 6, 3, 4, 1, 2, 
	8, 7, 6, 5, 4, 3, 2, 1
};


/*
 * local stuff
 */

#ifdef WITH_MOLCAS
int molcas_init_hook(int flags);
#endif

#ifdef WITH_GUK
int guk_init_hook(int flags);
#endif

#ifdef WITH_STUB
int stub_init_hook(int flags);
#endif

static int (*init_hooks[MAX_HOOKS]) (int flags);

int guk_phis_init(int *arg)
{
        char buffer[MAXLENGTH];
        int rc = 0;
	int module = (*arg & 0xFFFF);
	int flags = (*arg >> 16);

	if (flags & FLAG_RSVD1) {
		fprintf(stderr,
			"The functionality to display frozen/deleted orbitals\n"
			"has been removed from the PHIS library.\n"
			"The authors feel that if you want information about\n"
			"certain orbitals you should not discard them during\n"
			"the transformation.\n"
			"However, if you think this functionality is absolutely\n"
			"necessary for your application you can convince me to\n"
			"put it back in: Victor Vysotskiy <victor.vysotskiy@pci.uni-heidelberg.de>\n"
			"Terminating.\n");
		exit(0);
	}

#ifdef WITH_MOLCAS
        init_hooks[MOLCAS] = &molcas_init_hook;
#endif

#ifdef WITH_GUK
        init_hooks[GUK] = &guk_init_hook;
#endif

#ifdef WITH_STUB
        init_hooks[STUB] = &stub_init_hook;
#endif

	if (module >= MAX_HOOKS || init_hooks[module] == NULL)
		return(rc);

	/* call the init_hook */
	rc = init_hooks[module] (flags);
	return(rc);
}

void guk_phis_get_info(int *nSym, int *nBas, int *nAtoms)
{
	guk_interface.info(nSym, nBas, nAtoms);
}

void guk_phis_get_epsi(double *E_hf, double *e, int *n)
{
  guk_interface.epsi(E_hf, e, n);
}

void guk_phis_get_sym(int *s, int *n)
{
	guk_interface.sym(s, n);
}

void guk_phis_get_occ(double *o, int *n)
{
	guk_interface.occ(o, n);
}

void guk_phis_get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint)
{
	guk_interface.next_Vpqrs(p, q, r, s, vint);
}

void guk_phis_get_scfvec(double *C, int *n, int *len)
{
	guk_interface.scfvec(C, n, len);
}

void guk_phis_get_overlap(double *S, int *n)
{
	guk_interface.overlap(S, n);
}

void guk_phis_get_geometry(int *n, double *geo, double *Z_nuc, double *E_nuc)
{
	guk_interface.geometry(n, geo, Z_nuc, E_nuc);
}


void guk_phis_get_ao(int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
		 int *nmb_ao,int *max_nmb_cc)
{
	guk_interface.ao(polynomial,nmb_cc,cc,alpha,center,
		     nmb_ao,max_nmb_cc);
}

void guk_phis_list_active(int *list, int *n)
{
	guk_interface.list_active(list, n);
}

void guk_phis_get_dip(double *x, double *y, double *z, int *n)
{
	guk_interface.dip(x, y, z, n);
}

void guk_phis_get_vel(double *x, double *y, double *z, int *n)
{
	guk_interface.vel(x, y, z, n);
}

void guk_phis_get_quad(double *xx, double *xy, double *xz,
                   double *yy, double *yz, double *zz, int *n)
{
	guk_interface.quad(xx, xy, xz, yy, yz, zz, n);
}

void guk_phis_get_oneel(double *h, double *t, int *n)
{
	guk_interface.oneel(h, t, n);
}

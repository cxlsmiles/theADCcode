#include <stdlib.h>
#include <stdio.h>

#include "molcas.h"
#define MAXWORK  2048

int 
molcas_init_hook(int flags)
{
  printf("inti molcas\n");
  /* FIXME: allocate this locally in m_setup(). */
        int len = MAXWORK;
        double work[MAXWORK];

	interface.info = &phis_mc_info;
	interface.epsi = &phis_mc_epsi;
	interface.sym = &phis_mc_sym;
	interface.occ = &phis_mc_occ;
	interface.next_Vpqrs = &phis_mc_nextint;
	interface.scfvec = &phis_mc_scfvec;
	interface.overlap = &phis_mc_overlap;
	interface.ao = &phis_mc_ao;
	interface.geometry = &phis_mc_geometry;
	interface.list_active = &phis_mc_active;
	interface.dip = &phis_mc_dip;
	interface.vel = &phis_mc_vel;
	interface.quad = &phis_mc_quad;
	interface.oneel = &phis_mc_oneel;

	if (flags & SYM_BLOCKED) {
		interface.cap =
			(HAVE_INFO | HAVE_EPSI | HAVE_SYM | HAVE_OCC | 
			 HAVE_INT | HAVE_NEXTINT |
			 HAVE_AO | HAVE_GEOMETRY | HAVE_LOA  |
			 HAVE_OVERLAP | HAVE_DIP | HAVE_QUAD | HAVE_ONEEL ) ;
	} else {
		interface.cap = 
			(HAVE_INFO | HAVE_EPSI | HAVE_SYM | HAVE_OCC |
			 HAVE_GEOMETRY | 
			 HAVE_OVERLAP | HAVE_DIP | HAVE_QUAD | HAVE_ONEEL );
	}

	if (mc_query_vel() == 0)
		interface.cap |= HAVE_VEL;

	/* initialize some data */
	mc_setup(&flags, work, &len);

	return interface.cap;
}

/* The FORTRAN Stop directive returns with EXIT_SUCESS which is not
   what I want. */
void
mc_stop(void)
{
	exit(EXIT_FAILURE);
}

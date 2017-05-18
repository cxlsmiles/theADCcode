#include <stdio.h>
#include <stdlib.h>
#include "phis.h"
#include "config.h"

int
main(int argc, char **argv) {

	int p,q,r,s;
	double Vpqrs;
	int c = 0; 
#ifdef  WITH_GUK
	int backend = GUK;
#endif
#ifdef  WITH_MOLCAS
        int backend = MOLCAS |(SYM_BLOCKED << 16);
#endif        
        long cap;
        cap = phis_init(&backend);
	if ( !(cap & HAVE_NEXTINT) ) {
		fprintf(stderr, "No integrals!\n");
		exit(-1);
	}
      printf("   p   q   r   s |    Vpqrs\n"
	     "-----------------+-------------\n");

      phis_get_next_Vpqrs(&p,&q,&r,&s,&Vpqrs);
      do {
	      if (p==s && q==r && Vpqrs<0) {
		      printf("Problem with integral (%i,%i|%i,%i): %f\n",
			     p,q,r,s, Vpqrs);
		      exit(1);
	      }
	      if (p==s && q==r && p==q ) {
	printf(" %3d %3d %3d %3d | %12.8f\n",p,q,r,s,Vpqrs);
	      }
	      phis_get_next_Vpqrs(&p,&q,&r,&s,&Vpqrs);
      } while (s > 0);

      return 0;
}


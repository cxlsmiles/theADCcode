#include "./guk.h"

static char *rcsid = "$Id: phis_guk_nextint.c,v 1.2 2007/06/28 15:58:55 redo Exp $";

/* makes the data in Guk_buffer available as a mainfile block structure */
guk_mainfile_block *Buf = (guk_mainfile_block *)Guk_buffer;

/* makes the data in Guk_buffer available as a mainfile block structure (redo for huge basis set)*/
guk_mainfile_block_huge *Buf_huge = (guk_mainfile_block_huge *)Guk_buffer;

/* these preserve the state across calls */
static FILE *main_fp;
static int x = 0;

void
phis_guk_next_Vpqrs(int *p, int *q, int *r, int *s, double *Vpqrs) 
{

	if (! main_fp) { /* first call */
		if( (main_fp = fopen(mainfile_name, "rb")) == NULL ) {
			perror("Error opening mainfile"); 
			exit (1); 
		}
	}

	if (x == 0) { /* refill buffer */
		if (! guk_get_block(main_fp, 0)) {
			fclose(main_fp);
			*s = 0; main_fp = NULL;
			return;
		}
	}

	/* we want p>q and r>s */
        // redo for big basis set
        *q = 0;
        *p = 0;
        *s = 0;
        *r = 0;

        if (nBas > 256) {
          // for huge basis set
	  memcpy (q, Buf_huge->index[x].i, 2);
          memcpy (p, Buf_huge->index[x].j, 2);
          memcpy (s, Buf_huge->index[x].k, 2);
          memcpy (r, Buf_huge->index[x].l, 2);

          *Vpqrs = Buf_huge->integral[x];
	  
	  /* buffer exhausted ? */
	  if (++x == Buf_huge->integral_count) x = 0;
        }
        else {
          // for standard basis set
	  *q = Buf->index[x].i;
          *p = Buf->index[x].j;
          *s = Buf->index[x].k;
          *r = Buf->index[x].l;
  	  
          *Vpqrs = Buf->integral[x];

	  /* buffer exhausted ? */
	  if (++x == Buf->integral_count) x = 0;
        }

	return;
}

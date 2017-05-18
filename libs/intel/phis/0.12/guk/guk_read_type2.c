#include "phis.h"
#include "guk.h"

#include <math.h>
#include <string.h>
#include <assert.h>

/* For GUK version 5
 * every block in a section of type 2 has the following structure:
 *
 *   double nucl_eng;                        nuclear energy
 *   double dx,dy,dz;                        gradients
 *   double s[n],t[n],h[n],x[n],y[n],z[n];   parts of the matrices in packed storage
 *   ...
 *   int n;                                  nmb. of matrixele. in this blk.
 *   int i3;                                 ??? (always  3)
 *   int i48;                                ??? (always 48)
 *   int blk_no;                             block number
 *   int blk_len;                            block length
 *
 * In contrast to that, in GUK version 6 the matrices are stored consecutively.
 * Every matrix starts at a block boundary.
 *
 *
 * input:
 *   double *buf      where to store the result
 *   int buflen       the length of the result buffer
 *   int what         see enum type2 in guk.h for valid values
 */

void 
guk5_read_type2(double *buf, int buflen, int what)
{
	int i, n, ind;
	FILE *fp;
	int sec, pos, len;

	int wc;

	double *ptr = (double *) Guk_buffer;

	if ((fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		perror( "Cannot open dumpfile" ); 
		exit(EXIT_FAILURE);
	}

	sec = guk_sec_by_type(CORE_AO);
	pos = Summary.entry[sec].pos;
	len = Summary.entry[sec].len;

	switch (what) {
	case ENUC :
		wc = guk_get_block(fp, pos);
		buf[0] = ptr[0];
		break;
	case OLAP :
	case KIN :
	case CORE :
	case DIPX :
	case DIPY :
	case DIPZ :
		ind = 0;
		for(i=0; i<len; i++) {
			(void) guk_get_block(fp,pos+i);
			n = ((int*)Guk_buffer)[1019];
			memcpy(buf+ind, ptr+4+n*what, (size_t)8*n);
			ind += n;
		}
		assert(ind==buflen);
		break;
	default :
		fprintf(stderr,
			"unknown matrix type: %i\n", what);
		exit(EXIT_FAILURE);
	}

	fclose(fp);
	return;
}

void
guk6_read_type2(double *buf, int buflen, int what)
{
	FILE *fp;
	int sec, pos;
	int wc;

	if ((fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		perror( "Cannot open dumpfile" ); 
		exit(EXIT_FAILURE);
	}

	sec = guk_sec_by_type(CORE_AO);
	pos = Summary.entry[sec].pos;

	if (what == ENUC) {
		double *ptr = (double *) Guk_buffer;

		wc = guk_get_block(fp, pos);
		buf[0] = ptr[0];
		return;
	} else 
		pos += 1;

	pos += what*(buflen/511+1);
	wc = guk_read_data(fp, (char *) buf, pos, 0);
	assert(wc == buflen);

	fclose(fp);
	return;
}

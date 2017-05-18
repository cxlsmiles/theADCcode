/*
 * Very simple program that dumps the content of the first block of a
 * given section. It assumes 4-byte integers!
 *
 */

#include "../guk.h"

#include <unistd.h>

#define NWORDS 511
#define OPTIONS "hs:p:"

void
print_usage(char *progname)
{
	printf("Usage: %s -s<num>\n"
	       "Dumps the first block of section <num>\n", progname);
	exit(0);
}


int
main(int argc, char **argv)
{
	int opt;

	char *dumpfile_name;
	FILE *dumpfile_fp;

	int sec, pos = 0;
	char buf[4096];

	int i, j, wc;
	int *intptr;
	double *dblptr;
	char *charptr;

	if (argc < 2) print_usage(argv[0]);
	while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
		switch (opt) {
		case 's' :
			sec = atoi(optarg);
			break;
		case 'p' :
			pos = atoi(optarg);
			break;
		case 'h' :
		default:
			print_usage(argv[0]);
		}
	}

	if ((dumpfile_name = getenv("ed3")) == NULL) {
		fprintf(stderr, "You must assign the name of the dumpfile to ed3.\n");
		exit( 1 );
	}
	if ((dumpfile_fp = guk_open_dump_file(dumpfile_name)) == NULL) {
		fprintf(stderr, "Cannot open dumpfile %s\n", dumpfile_name); 
		exit( 1 );
	}
	
	/* FIXME: be consisten about where the numbering starts for sec and pos */
	if (! pos) 
		pos = (int) Summary.entry[sec-1].pos;
	wc = guk_read_data(dumpfile_fp, buf, pos, NWORDS);

	printf("Dumping section %i, found at block %i\n", sec, pos);

	dblptr = (double *) buf;
	intptr = (int *) buf;
	charptr = (char *) buf;
	for (i=0; i<NWORDS; i++) {

		printf("%3i:\t", i);
		printf("% 12i  ", *intptr++);
		printf("% 12i\t", *intptr++);
		printf("% 12g\t", *dblptr);
		printf("%.8s\t", charptr);
		printf("%08X", *dblptr++);
		printf("\n");
		charptr += 8;
	}

	return(0);
}

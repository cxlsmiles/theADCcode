#include <string.h>

#include "./guk.h"

/************************** global variables *****************************/

#ifdef DEBUG
static int  Guk_debug = 10;         /* Verbose output if nonzero             */
#else
static int  Guk_debug = 0;
#endif
static int  Guk_wc;		   /* Word count of the last block read     */
static int  Guk_blkno;		   /* Block number of the last block        */

char Guk_buffer[4096];             /* global read buffer - holds one block 
				      at a time */
Type_0  Summary;                   /* set by guk_open_file */

FILE *
guk_open_dump_file(char *filename) 
/* Opens the GUK dumpfile specified by filename and tries to 
   read its summary.
   Returns the filepointer or NULL on error */
{
	int sec;
	FILE *fp;

	if ((fp = fopen(filename, "rb")) == NULL)
		return NULL;

	if ((fread(&Summary, sizeof(Summary), 1, fp)) != 1)
		return NULL;

	if (Guk_debug) {
		printf("Contents of Block Zero:\n");
		printf("Maxb = %i kblkla = %i WordCount = %i\n", 
			Summary.max_blocks, Summary.used_blocks, Summary.word_count);

		for (sec=0; sec<507; ++sec) {
			if (Summary.entry[sec].type)
				printf("Section: %3d\t Pos: %5d\t Type: %4d\t Len: %4d\n", 
				       sec+1, 
				       Summary.entry[sec].pos,
				       Summary.entry[sec].type,
				       Summary.entry[sec].len); 
		}
	}
	return fp;
}


void *
guk_get_section(FILE *fp, int type) 
{
	char *section;
	int word_count;
	int sec, len, pos;
	int i, size;
	int words_wanted, words_stored;

	/* look up the section */
	if ((sec = guk_sec_by_type(type)) == -1) return NULL;
	len = Summary.entry[sec].len;
	pos = (int) Summary.entry[sec].pos;
  
	/* allocate using init_section */
	section = guk_init_section(NULL,type);

	/* the size is stored after the data pointer in every section */
	size = (int) *((char**)section+1);

	/* round upwards next word => +7 */
	words_wanted = (size+7)/8; 

	/* this positions the file just in front of the first block in the section */
	guk_position_file(fp, pos);

	words_stored = 0;

	for( i=0; i<len; ++i ) {

		word_count = guk_get_block(fp, 0); /* read the next block */

		/* check if there is enough space in allocated memory */
		if(words_stored + word_count > words_wanted ) {
		  fprintf(stderr,
			  "guk_get_section: not enough space in Type_%d allocated for stored data\n"
			  "                 %d words needed > %d words (%d bytes) allocated\n",
			  type,words_stored+word_count,words_wanted,size);
		  exit(-1);
		}

		/* copy from read buffer to section->data buffer */
		memcpy(*((char**)section)+8*words_stored, Guk_buffer, 8*word_count);
		words_stored += word_count;

		/* break if Type structure filled */
		if(words_stored==words_wanted)
			break;
	}

	if (Guk_debug) {
		printf("Section of type %i requires %i words.\n", type, words_wanted);
		printf("no. words actually stored: %i\n", words_stored); 
	}

	/* need to initialize again for the assignment of the non-pointer vars */
	guk_init_section(section,type);

	return section;
}

int 
guk_read_data(FILE *fp, char *buf, int blkno, int n_words) 
/* reads n_words 8-byte words into buf starting at block number blkno or 
   at the next block if blkno=0
   reads across block boundaries                                            
   if n_words=0, reads as many consecutive words as it can find according 
   to the word count of the blocks                                
   Returns the number of words read                                         
   if n_words given, does not check the word count */
{
	int i;
	int n_blocks, left_over;
	int word_count;
	int words_stored = 0;
  
	if (blkno) guk_position_file(fp, blkno);

	switch (n_words) {

	case 0 : 
		while (1) { 
			word_count = guk_get_block(fp, 0); 

			if (Guk_debug)
				printf("word count: %i\n", word_count);

			/* copy from read buffer */
			memcpy(buf+8*words_stored, Guk_buffer, 8*word_count);
			words_stored += word_count;

			if (word_count < 511) break;
		}
		break;

	default :
		if (Guk_debug)
			printf("Reading %i consecutive words, ignoring word count\n", n_words);

		n_blocks = (int) (n_words/511);
		left_over = n_words - n_blocks*511;
		for (i=0; i<n_blocks; ++i) {
			guk_get_block(fp, 0); 
			word_count = 511;

			/* copy from read buffer */
			memcpy(buf+8*words_stored, Guk_buffer, 8*word_count);
			words_stored += word_count; 
		}

		/* read what's left over */
		guk_get_block(fp, 0);
		word_count = left_over;

		memcpy(buf+8*words_stored, Guk_buffer, 8*word_count);
		words_stored += word_count;
		break; 
	}

	return words_stored;
}

void
guk_getpara(FILE *dump_fp)
{
	/* the constants are determined using the size of the subsections 
	   (the continuous parts) of the sections type 1 and 3.
	   
	   subsection sizes of Type_1:
	     mxprim
	     mxprim*maxaqm + maxat
	     (80 + maxat*8)/8
	     (7 * mxshel + 2)/2

	   subsection sizes of Type_3:
	     (job specification ignored)
	     2*maxorb + 4
	*/

	int sec,pos;

	/* first use section type 1 */
	if ((sec = guk_sec_by_type(1)) < 0) {
		fprintf(stderr, "Cannot find section type 1 in dfile.\n");
		exit(1);
	}
	pos = Summary.entry[sec].pos;
	guk_position_file(dump_fp, pos);

	guk_mxprim = guk_sizeof_subsection(dump_fp);
	guk_maxaqm = guk_sizeof_subsection(dump_fp);
	guk_maxat  = guk_sizeof_subsection(dump_fp) - 10;
	guk_mxshel = (guk_sizeof_subsection(dump_fp) - 2)/7 * 2; /* FIXME: just works for 4-byte integers */
	
	/* calc correct maxaqm using mxprim and maxat */
	guk_maxaqm = (guk_maxaqm - guk_maxat)/guk_mxprim;
	

	/* now for section type 3 */
	if ((sec = guk_sec_by_type(3)) < 0) {
		fprintf(stderr, "Cannot find section type 3 in dfile.\n");
		exit(1);
	}
	pos = Summary.entry[sec].pos;
	guk_position_file(dump_fp, pos);

	guk_sizeof_subsection(dump_fp);
	guk_maxorb = (guk_sizeof_subsection(dump_fp) - 4)/2;

	if(Guk_debug) {
		printf("MXPRIM = %d\n",guk_mxprim);
		printf("MAXAQM = %d\n",guk_maxaqm);
		printf("MAXAT  = %d\n",guk_maxat);
		printf("MXSHEL = %d\n",guk_mxshel);
		printf("MAXORB = %d\n",guk_maxorb);
	}
}

/* ------------------- low level routines -------------------------------- */

void 
guk_position_file(FILE *fp, int pos) 
{
	if (fseek(fp, pos*BLOCKSIZE, SEEK_SET)) {
		perror("fseek");
		exit(1); 
	}
	return;
}

int 
guk_get_block( FILE *fp, int blkno )
/* low level read routine:                                                  */
/* read block with number blkno, or the next block if blkno = 0             */
/* Returns the word count                                                   */
{  
	unsigned long long last_word;

	if (blkno) guk_position_file(fp, blkno);

	if ((fread(Guk_buffer, sizeof(char), 4096, fp)) != 4096) {
		perror("fread");
		exit(1); 
	}

	/* this is tricky: (buf+OFFSET) is the starting address of the last word 
	   in the block, (unsigned long long *) makes it the starting address of 
	   an 64 bit integer and the final asterisk dereferences the pointer, 
	   returning an integer value */
	last_word = *(unsigned long long *) (Guk_buffer + 4088);

	Guk_wc = (int) last_word;
	last_word = last_word >> 32;
	
	/* Fixes block number: the number stored in the dumpfile blocks is 
	   always one unit greater than the position marked in the summary. */
	Guk_blkno = (int) last_word - 1; 

	if (Guk_debug)
		printf("Read block no. %i, word count = %i\n", Guk_blkno, Guk_wc);

	return Guk_wc;
}

int
guk_sec_by_type(int type) 
{
	int sec;
  
	if (Guk_debug)
		printf ("Looking up type %i\n", type);
	for (sec=0; sec<508; ++sec)
		if (Summary.entry[sec].type == type)
			break;

	/* section not found: this should also work in the case that the */
	/* desired type is in the last section                           */
	if (Summary.entry[sec].type != type) {
		printf("Section of type %i not found in dumpfile.\n", type);
		return -1;
	}

	if (Guk_debug)
		printf("Found in section %i\n", sec+1);

	return sec;
}

int
guk_sizeof_subsection(FILE *fp)
{
	int word_count,words_stored = 0;

	do {
		word_count = guk_get_block(fp,0);
		words_stored += word_count;
	} while(word_count==511);

	return(words_stored);
}

/**********************	other helper functions *************************/

void
copy_and_expand(double *M, double *buf, int n)
{
/* 
 * Takes a matrix in packed storage mode in buf 
 * and expands it into a full n*n matrix in M.
 */
	int idx = 0;
	int i,j;
	double x;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			x = buf[idx++];
			M[i*n+j] = x;
			M[j*n+i] = x;
		}
	}

	return;
}

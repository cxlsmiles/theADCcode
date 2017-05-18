/* $Id: typedefs.h,v 1.2 2007/06/28 15:58:55 redo Exp $ 
 *
 * Internal data structures for the GUK dumpfile
 *
 * a word is 8 bytes
 * a block is 512 words (i.e. 4096 bytes)
 * a section streches over one or more blocks  
 *
 */

/* GUK section types */
#define SUMMARY               0
#define SCF_DATA              3
#define BASIS_SET_DATA        1
#define GEOMETRY_DATA        15
#define SYMMETRY_DATA        51
#define CORE_AO               2
#define CORE_MO            1004

// REDO FOR GUK BIG BASIS

#define MAXINT      340  /* maximum number of integrals per block */
#define MAXINTHUGE  255 
#define BLOCKSIZE  4096

/*** Type_0: the summary section *******************************************/
typedef struct {
  struct { 
    short unused;		/* not used */
    short pos;		        /* block BEFORE this section */
    short type;		        /* type of section */
    short len;		        /* length in blocks */
  } entry[508];			/* max. 509 sections may be used */
  int          max_blocks;	/* maximum number of blocks */
  int          used_blocks;	/* actual number of blocks used */
  char unused[2][8];		/* two unused words */
  int block_number;
  int word_count;          	/* number of words used, always 509 */ 
} Type_0;


/*** Type_1: the basis set ************************************************/
typedef struct {
  char *data;         /* here all the data is stored, the rest points to this part */
  int size;           /* size of data segment */

  double *ex;         /* #mxprim          list of unique primitive exponents */
  double *cc;         /* #(maxaqm*mxprim) list of prim. contraction coefficients */
  double *z;          /* #maxat           atomic number of all atoms */

  char *job_title;    /* #80              job description given by user */
  char *atom_labels;  /* #(maxat*8)       right aligned atomic labels */

  int *kstart;        /* #mxshel          definition of AOs */
  int *katom;         /* #mxshel          (see phis_guk_ao for description) */
  int *ktype;         /* #mxshel */
  int *kng;           /* #mxshel */
  int *kloc;          /* #mxshel */
  int *kmin;          /* #mxshel */
  int *kmax;          /* #mxshel */

  int n_shells;       /* total number of shells */
  int n_atoms;        /* total number of atoms */
  int n_basis;        /* number of basis functions */
} Type_1; 


/*** Type_3: the vectors section *******************************************/
typedef struct {
  char *data;         /* here all the data is stored, the rest points to this part */
  int size;           /* size of data segment */

  char *user;         /* #8    the job specifications */
  char *date;         /* #8    */
  char *time;         /* #8    */
  char *prog;         /* #8    */
  char *scftype;      /* #8    */
  char *acct;         /* #8    */
  char *more;         /* #104  */
  char *job_title;    /* #80   */          

  double *orb_eng;    /* #maxorb  orbital energies */    
  double *occ_num;    /* #maxorb  occupation numbers */
  double ehf_tot;     /*          total hartree fock energy */
  int n_gtos;         /*          number of gaussian type orbitals */
  int nlcbf;          /*          ??? */
  int n_vec;          /*          number of scf-vectors */
  int *unused;        /* #3       ??? */

  int *ilifc;          /* #maxorb    ??? */
  int *ntran;          /* #maxorb    number of sym. equivalent AOs for this scfvec */
  int *itran;          /* #maxorb*3  indices of these AOs */
  double *ctran;       /* #maxorb*3  and coefficients for the linear combination */
  int iftran;          /*            ??? */
} Type_3; 
/* variable number of blocks containing the SCF vectors follow */


/*** Type_51: the symmetry information ************************************/
typedef struct {
  char *data;          /* here all the data is stored, the rest points to this part */
  int size;            /* size of data segment */

  int n_irrep;         /*          order of the point group */
  int *mult_table;     /* #8*8     multiplication table (actually constant) */
  int *symlabels_ao;   /* #maxorb  symmetry lables in the AO basis */
  int *symlabels_mo;   /* #maxorb  same in the MO basis */
} Type_51; 


/*** Type_2: core properties in the AO basis *****************************
 *
 * Actually I do not use a structure for this, as the data seems to be
 * rather block-oriented. I THINK the structure of this section is like
 *
 * 1st block: double potnuc, double pot_grad[3]
 * 2nd block: ovelap
 * 3rd block: kinetic energy
 * 4th block: core Hamiltonian
 * 5th block: x coordinates of the nuclei
 * 6th block: y coordinates
 * 7th block: z coordinates
 *
 * This, however, is only a working hypothesis; and I have no idea what happens
 * if you use more than 23 basis functions...
 *
 */

/*** Type_1005: contains the core/active lists *************************
 *
 * Note that this section is overwritten on a restart!
 */
typedef struct {
  char *data;          /* here all the data is stored, the rest points to this part */
  int size;            /* size of data segment */

  int nAct;
  int *list_of_act;    /* #maxorb */
  int act_print_flag;
  int nFro;
  int *list_of_fro;    /* #maxorb */
  int fro_print_flag;
/* Here follows a lot of junk giving this section a total length of 9
   blocks. However, that data seems to be uninitialized! 
   (cf. master.m, subroutine start2, around line 9131) */
} Type_1005;


/*** format of the integrals in the mainfile ****************************/
typedef struct {
  double integral[MAXINTHUGE];
  // REDO guk FOR BIG BASIS
  struct {
    unsigned char i[2];
    unsigned char j[2];
    unsigned char k[2];
    unsigned char l[2];
  } index[MAXINTHUGE];
  int integral_count;		/* Nr of integrals in the block */
  int unknown;		        /* ??? */
  int block_count;	
  int word_count;
} guk_mainfile_block_huge;

typedef struct {
  double integral[MAXINT];
  // REDO guk FOR BIG BASIS
  struct {
    unsigned char i;
    unsigned char j;
    unsigned char k;
    unsigned char l;
  } index[MAXINT];
  int integral_count;		/* Nr of integrals in the block */
  int unknown;		        /* ??? */
  int block_count;	
  int word_count;
} guk_mainfile_block;

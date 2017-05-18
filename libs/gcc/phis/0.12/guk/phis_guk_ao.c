#include "./guk.h"

#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795029L

/* phis_guk_ao:
 * reads atomic orbital information
 * output is close to the old "Hondo" format (see below)
 *
 * input:
 * ------
 * nmb_ao:     number of atomic orbitals expected
 * max_nmb_cc: maximal number of contractions
 *
 *
 * output:
 * -------
 * polynomial: (nmb_ao x 3)-matrix of the polynomial exponents
 * nmb_cc:     nmb_ao-list with number of contractions 
 * cc:         (nmb_ao x max_nmb_cc)-matrix of contraction coefficients
 * alpha:      (nmb_ao x max_nmb_cc)-matrix of alpha's
 * center:     nmb_ao-list of centers
 *
 *
 * information in section type 1:
 * ------------------------------
 * ex        = list of unique primitive exponents (identical functions
 *             on equivalent centers not repeated). See below.
 * cc        = list of unique contraction coefficients (identical
 *             functions on equivalent centers not repeated). See below.
 *             The coefficients are in terms of normalized primitives.
 * z         = atomic number of all atoms (including equivalent).
 * job_title = title of last scf run.
 * atom_labels = label (right aligned) of all atoms (including equivalent).
 * kstart(i) = location of first primitive of shell i (i=1,nshell)
 *             in the list of unique exponents (ex) and contraction
 *             coefficients (cc) The contraction coefficients for shell
 *             i are in cc(j,k), j=kstart(i),kstart(i)+kng(i)-1,
 *             where k numbers the actual l+1 values for the shell.
 *             The exponents are in ex(j).
 * katom(i)  = center where shell i resides, i=1,nshell.
 * ktype(i)  = maximum l+1 value for shell i (1=s, 2=p, 3=d, 4=f).
 * kng(i)    = # of primitives in shell i, i=1,nshell
 * kloc(i)   = index of first b.f. of shell i in the complete list of
 *             basis functions. The bf's for this shell go from
 *             kloc(i) to kloc(i)+kmax(i)-kmin(i)+1 (see below).
 * kmin(i),  = basis functions in shell i are of type kmin(i) to kmax(i).
 * kmax(i)     Types are:
 *             s=1,
 *             px=2, py=3, pz=4,
 *             dxx=5, dyy=6, dzz=7, dxy=8, dxz=9, dyz=10,
 *             fxxx=11, fyyy=12, fzzz=13, fxxy=14, fxxz=15, fxyy=16.
 *             fyyz=17, fxzz=18, fyzz=19, fxyz=20.
 * n_shells  = total number of shells. A shell is a contraction or group
 *             of contractions with identical exponential part, e.g. a p
 *             shell (three basis functions), a split valence sp shell
 *             (four basis functions), etc. This counts the total, i.e.
 *             including identical shells on equivalent centers.
 * n_atoms   = total number of atoms (including equivalent).
 * n_basis   = total number of basis functions: s+3*p+4*sp+6*d+10*f.
 *
 */
void
phis_guk_ao(int *polynomial,int *nmb_cc,double *cc,double *alpha,int *center,
	    int *nmb_ao,int *max_nmb_cc)
{
  FILE *fp;
  Type_1 *section1;
  int thisshell,typeindex,l_value;
  int i,j;

  /* conversion from type-index to exponents of x,y and z */
  int l_tab[3][20]={{0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,1,0,1,0,1},
		    {0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1},
		    {0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,0,1,2,2,1}};

  /* f_tab[i]: prefactor from angular momentum for normalization */
  double f_tab[20]={ 1.,
		     1./2,  1./2,  1./2, 
		     3./4,  3./4,  3./4,
		     1./4,  1./4,  1./4,
		    15./8, 15./8, 15./8,
		     3./8,  3./8,  3./8, 3./8, 3./8, 3./8,
		     1./8};


  /* check, if there is enough space for all basis functions */
  if(*nmb_ao < section3->n_gtos) {
    *nmb_ao = -section3->n_gtos;
  }
  
  /* now read section of type 1 (basisset information) */
  if((fp=guk_open_dump_file(dumpfile_name))==NULL) {
    fprintf(stderr, "Cannot open dumpfile %s\n", dumpfile_name); 
    exit( 1 );
  }
  section1 = guk_get_section(fp,1);
  fclose(fp);

  /* check if there is enough space for the exponents */
  for(i=0;i<section1->n_shells;i++) 
    if(section1->kng[i] > abs(*max_nmb_cc)) 
      *max_nmb_cc = -section1->kng[i];

  /* if out of space: return space requirements as negative values */
  if(*max_nmb_cc < 0 || *nmb_ao < 0) 
    return;


  /* forall basis functions do */
  for(i=0;i<section3->n_gtos;i++) {
    
    /* at first we need the shell of the basis function */
    for(j=0;j<section1->n_shells;j++) {

      if(section1->kloc[j]-1<=i && 
	 i<=section1->kloc[j]-1+section1->kmax[j]-section1->kmin[j]) {

	thisshell=j;
	break;
      }
    }

    /* get the type-index: s=0,px=1,py=2 ... */
    typeindex = i - (section1->kloc[thisshell]-1) + section1->kmin[thisshell]-1;

    /* convert type-index to polynomial */
    l_value=0;
    for(j=0;j<3;j++) {
      polynomial[3*i+j] = l_tab[j][typeindex];
      l_value += l_tab[j][typeindex];
    }
    
    /* get number of primitives for this basisfunction */
    nmb_cc[i] = section1->kng[thisshell];

    for(j=0;j<nmb_cc[i];j++) {
      /* get exponents */
      alpha[*max_nmb_cc*i+j] = section1->ex[section1->kstart[thisshell]-1+j];

      /* and normalized contraction coefficients */
      cc[*max_nmb_cc*i+j]
	= section1->cc[l_value*guk_mxprim + section1->kstart[thisshell]-1+j]
	/(pow(PI,3./4)*sqrt(f_tab[typeindex])
	  /pow(2*alpha[*max_nmb_cc*i+j],(3./2+l_value)/2));
    }


    /* at least the number of the center where shell i resides */
    center[i] = section1->katom[thisshell];
  }
}

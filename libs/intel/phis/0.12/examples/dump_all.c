
/* dump_all.c
 * 
 * Small sample application, which dumps all the available information
 * of the SCF calculation to stdout.
 *
 * The same thing in FORTRAN can be found in 'fdump_all.f'.
 *
 * The different interfaces are specified using flags (see usage). The
 * SCF files are specified using environmental variables (e.g. ed3 for
 * the guk dumpfile)
 *
 * After some startup checks the following sections of this program are
 * always structured like
 *
 *   1.) Check if section is supported using the variable 'cap' and 
 *       'HAVE_<section>'
 *   2.) The needed space is allocated
 *   3.) The library function is called
 *   4.) Error checking is done:
 *         Does the information fit into the arrays?
 *         (Since the arrays are allocated using the information of the 
 *         'phis_get_info' function this should always be the case)
 *   5.) (formatted) printing of the info to stdout
 *   6.) Clean up (arrays and vectors are freed)
 *
 * (This structure is pointed out in the AO section, see below)
 *
 * JB(3/02)
 */

#include "phis.h"
#include <math.h>

#define OPTIONS "mgis"

void print_matrix(double * matrix,int nmb_rows,int nmb_cols);


/*----------------------------------------------------------------
 * usage
 */
void print_usage(char *progname)
{
  printf("%s [-m|-g] [-ias] dumps available information\n"
	 "-m\t selects MOLCAS backend\n"
	 "-g\t selects GAMESS-UK backend\n"
	 "-i\t prints two electron integrals (handle with care)\n"
	 "-s\t block by symmetry\n",
	 progname);

  exit(-1);
}


int main(int argc,char **argv)
{

  long cap;
  int backend, flags = 0;
  int opt,print_integrals=0;
  int nSym, nBas, nAtoms, nAct;

  int n;

  int i,j;

  /*----------------------------------------------------------------
   * check out the arguments/options given
   */
  if (argc < 2) print_usage(argv[0]);
  while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
    switch (opt) {
    case 'm' :
      backend = MOLCAS;
      break;
    case 'g' :
      backend = GUK;
      break;
    case 'i' :
      print_integrals=1;
      break;
    case 's' :
	    flags |= SYM_BLOCKED;
	    break;
    default:
      print_usage(argv[0]);
    }
  }

  /*----------------------------------------------------------------
   * some start up info 
   */
  {
	  int maj, min, pat;
	  phis_version(&maj, &min, &pat);
	  printf("dump_all using PHIS version %i.%i.%i\n"
		 "------------------------------------\n"
		 "(prints out nearly all information, which is\n"
		 "implemented for the specified interface)\n",
		 maj, min, pat);
  }

  backend |= (flags << 16);
  cap = phis_init(&backend);

  /*----------------------------------------------------------------
   * print out the info 
   */
  if ( cap & HAVE_INFO ) {
    printf("\n<info>\n\n");

    phis_get_info(&nSym, &nBas, &nAtoms);
    
    printf("no. of irreps:          %i\n", nSym);
    printf("no. of basis functions: %i\n", nBas);
    printf("no. of atoms:           %i\n", nAtoms);

    printf("\n</info>\n\n");
  }
  else {
    fprintf(stderr,"%s needs at least the 'phis_get_info' function to work\n",
	    argv[0]);
    exit(-1);
  }

  /*----------------------------------------------------------------
   * print the list of active orbitals
   */
  if (cap & HAVE_LOA) {
	  int n = nBas;
	  int *list;

	  printf("<!--------------------------------------------------------------"
		 "-------------->\n\n"
		 "\n<active>\n\n");
	  list = (int *) malloc_we(n*sizeof(int));
	  phis_list_active(list, &n);

	  if (n<0) {
		  fprintf(stderr, "Error reading active list.\n");
		  exit(-1);
	  }

	  for (i=0; i<n; i++) printf(" %4d %4d\n", i+1, list[i]);

	  free(list);
	  printf("\n</active>\n\n");
  }

  /*----------------------------------------------------------------
   * print energies 
   */
  if ( cap & HAVE_EPSI ) {
    int n=nBas;
    double e_hf,*epsi;

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<epsi>\n\n");

    epsi=(double*)malloc_we(n*sizeof(double));
  
    phis_get_epsi(&e_hf,epsi,&n);

    if(n!=nBas) {
      fprintf(stderr,"Sorry, not enough space to save the energies\n");
      exit(-1);
    }

    printf("total hartree fock energy: %10.4f\n"
	   "orbital energies following:\n",e_hf);
    for(i=0;i<n;i++) printf(" %4d %8.4f\n",i+1,epsi[i]);

    free(epsi);

    printf("\n</epsi>\n\n");
  }

  /*----------------------------------------------------------------
   * print occupation numbers 
   */
  if ( cap & HAVE_OCC ) {
    int n=nBas;
    double *occ;

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<occ>\n\n");

    occ=(double*)malloc_we(n*sizeof(double));
  
    phis_get_occ(occ,&n);

    if(n!=nBas) {
      fprintf(stderr,
	      "Sorry, not enough space to save the occupation numbers\n");
      exit(-1);
    }

    printf("occupation numbers:\n");
    for(i=0;i<n;i++) printf(" %4d %8.4f\n",i+1,occ[i]);

    free(occ);

    printf("\n</occ>\n\n");
  }

  /*----------------------------------------------------------------
   * print symmetry information 
   */

  if ( cap & HAVE_SYM ) {
    int i,n=nBas;
    int *sym;

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<sym>\n\n");

    sym=(int*)malloc_we(n*sizeof(int));
  
    phis_get_sym(sym,&n);

    if(n!=nBas) {
      fprintf(stderr,
	      "Sorry, not enough space to save the symmetry labels\n");
      exit(-1);
    }

    printf("symmetry labels:\n");
    for(i=0;i<n;i++) printf(" %4d %4d\n",i+1,sym[i]);

    free(sym);

    printf("\n</sym>\n\n");
  }

  /*----------------------------------------------------------------
   * print geometry
   */
  
  if ( cap & HAVE_GEOMETRY ) {
    int n=nAtoms;
    int i,j;
    double *geom,*Z_nuc,E_nuc;

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<geometry>\n\n");

    geom =(double*)malloc_we(3*n*sizeof(double));
    Z_nuc=(double*)malloc_we(n*sizeof(double));

    phis_get_geometry(&nAtoms,geom,Z_nuc,&E_nuc);

    printf("Number of atoms: %d\n",n);
    printf("Nuclear energy : %g\n\n",E_nuc);

    printf(" charge |    x    |    y    |    z\n"
	   " -------+---------+---------+---------\n");
    for(i=0;i<n;i++) 
      printf(" % 5.2f  | % 6.4f   % 6.4f   % 6.4f\n",
	     Z_nuc[i],geom[3*i],geom[3*i+1],geom[3*i+2]);

    printf("\n</geometry>\n\n");
  }

  /*----------------------------------------------------------------
   * print AO information 
   */


  /* 1.) Check: information available? */
  if ( cap ) {

    /* some temporary variables are needed */
    int nmb_ao;
    int max_nmb_cc;

    int i;

    int *nmb_cc;
    int *polynomial;
    double *cc;
    double *alpha;
    int *center;
    
    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<ao>\n\n");

    /* 1a.) query dimensions */
    nmb_ao = 0; max_nmb_cc = 0;
    phis_get_ao(NULL,NULL,NULL,NULL,NULL,&nmb_ao,&max_nmb_cc);
    nmb_ao = abs(nmb_ao); max_nmb_cc = abs(max_nmb_cc);

    /* 2.) allocate some memory for the arrays */
    nmb_cc=(int*)malloc_we(nmb_ao*sizeof(int));
    polynomial=(int*)malloc_we(nmb_ao*3*sizeof(int));
    cc=(double*)malloc_we(nmb_ao*max_nmb_cc*sizeof(double));
    alpha=(double*)malloc_we(nmb_ao*max_nmb_cc*sizeof(double));
    center=(int*)malloc_we(nmb_ao*sizeof(int));

    /* 3.) ask the interface for the AOs */
    phis_get_ao(polynomial,nmb_cc,cc,alpha,center,&nmb_ao,&max_nmb_cc);
	
    /* 4.) everything ok? */
    if(nmb_ao<=0 || max_nmb_cc<=0) {
      fprintf(stderr,"error: phis_get_ao needs space for %d atomic orbitals\n"
	      "       and %d contraction coefficients\n",
	      abs(nmb_ao),abs(max_nmb_cc));
      exit(-1);
    }

    /* 5.) print AO table, the first three primitives are shown */
    printf("  AO |center| i j k | nmb_cc | first three alphas      |"
	   " first three cc's\n-----+------+-------+--------+-------"
	   "------------------+----------------------\n"); 
    for(i=0;i<nmb_ao;i++) {
      printf("%4d |%4d  |%2d%2d%2d |%5d   |%8.3f%8.3f%8.3f |%7.3f%7.3f%7.3f\n",
	     i+1,center[i],polynomial[3*i],polynomial[3*i+1],polynomial[3*i+2],
	     nmb_cc[i],
	     alpha[max_nmb_cc*i],
	     nmb_cc[i]>1 ? alpha[max_nmb_cc*i+1] : 0, /* print zeros if less */
	     nmb_cc[i]>2 ? alpha[max_nmb_cc*i+2] : 0, /* contractions */
	     cc[max_nmb_cc*i],
	     nmb_cc[i]>1 ? cc[max_nmb_cc*i+1] : 0,    /* print zeros if less */
	     nmb_cc[i]>2 ? cc[max_nmb_cc*i+2] : 0);   /* contractions */
    }

    /* 6.) free the used memory */
    free(nmb_cc);
    free(polynomial);
    free(cc);
    free(alpha);
    free(center);

    printf("\n</ao>\n\n");
  }
  
  /*----------------------------------------------------------------
   * print overlap matrix 
   */
  if ( cap & HAVE_OVERLAP ) {
    int n=0;
    double *s_matrix;

    phis_get_overlap(NULL,&n);
    n = abs(n);
    printf("%d\n",n);
    s_matrix=(double*)malloc_we(n*n*sizeof(double));
  
    phis_get_overlap(s_matrix,&n);

    if (n < 0) {
      fprintf(stderr,"Sorry, not enough space to save the S-matrix\n");
      exit(-1);
    }

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<overlap>\n\n");

    print_matrix(s_matrix,n,n);

    free(s_matrix);
    printf("\n</overlap>\n\n");
  }

  /*----------------------------------------------------------------
   * print scf vectors 
   */
  if ( cap & HAVE_SCFVEC ) {
    int nmb_scf, len_scf;
    int i,j,k;
    double *scfvectors;

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<scfvec>\n\n");

    /* query size of the SCF vectors array */
    nmb_scf = 0;
    len_scf = 0;
    phis_get_scfvec(NULL,&nmb_scf, &len_scf);

    /* allocate memory */
    nmb_scf = abs(nmb_scf);
    len_scf = abs(len_scf);
    scfvectors=(double*)malloc_we(nmb_scf*len_scf*sizeof(double));

    /* get vectors */
    phis_get_scfvec(scfvectors,&nmb_scf, &len_scf);

    /* and print them */
    print_matrix(scfvectors,len_scf,nmb_scf);

    free(scfvectors);
    printf("\n</scfvec>\n\n");
  }

  /*----------------------------------------------------------------
   * print dipole integrals (length form)
   */
  if ( cap & HAVE_DIP ) {
    int n=0;
    double *x, *y, *z;

    phis_get_dip(NULL, NULL, NULL, &n);
    n = abs(n);

    x = (double*)malloc_we(n*n*sizeof(double));
    y = (double*)malloc_we(n*n*sizeof(double));
    z = (double*)malloc_we(n*n*sizeof(double));
  
    phis_get_dip(x,y,z,&n);

    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<dipole>\n\n");

    printf("X component:\n");
    print_matrix(x,n,n);

    printf("Y component:\n");
    print_matrix(y,n,n);

    printf("Z component:\n");
    print_matrix(z,n,n);

    free(x); free(y); free(z);
    printf("\n</dipole>\n\n");
  }

  /*----------------------------------------------------------------
   * print dipole integrals (velocity form)
   */
  if ( cap & HAVE_VEL ) {
	  int n=0;
	  double *x, *y, *z;

	  phis_get_vel(NULL, NULL, NULL, &n);
	  n = abs(n);

	  x = (double*)malloc_we(n*n*sizeof(double));
	  y = (double*)malloc_we(n*n*sizeof(double));
	  z = (double*)malloc_we(n*n*sizeof(double));
  
	  phis_get_vel(x,y,z,&n);

	  printf("<!--------------------------------------------------------------"
		 "-------------->\n\n"
		 "<velocity>\n\n");

	  printf("X component:\n");
	  print_matrix(x,n,n);

	  printf("Y component:\n");
	  print_matrix(y,n,n);

	  printf("Z component:\n");
	  print_matrix(z,n,n);

	  free(x); free(y); free(z);
	  printf("\n</velocity>\n\n");
  }

  /*----------------------------------------------------------------
   * print quadrupole integrals
   */
  if ( cap & HAVE_QUAD ) {
	  int n=0;
	  double *xx, *xy, *xz, *yy, *yz, *zz;

	  phis_get_quad(NULL, NULL, NULL, NULL, NULL, NULL, &n);
	  n = abs(n);

	  xx = (double*)malloc_we(n*n*sizeof(double));
	  xy = (double*)malloc_we(n*n*sizeof(double));
	  xz = (double*)malloc_we(n*n*sizeof(double));
	  yy = (double*)malloc_we(n*n*sizeof(double));
	  yz = (double*)malloc_we(n*n*sizeof(double));
	  zz = (double*)malloc_we(n*n*sizeof(double));
  
	  phis_get_quad(xx, xy, xz, yy, yz, zz, &n);

	  printf("<!--------------------------------------------------------------"
		 "-------------->\n\n"
		 "<quadrupole>\n\n");

	  printf("XX component:\n");
	  print_matrix(xx,n,n);

	  printf("XY component:\n");
	  print_matrix(xy,n,n);

	  printf("XZ component:\n");
	  print_matrix(xz,n,n);

	  printf("YY component:\n");
	  print_matrix(yy,n,n);

	  printf("YZ component:\n");
	  print_matrix(yz,n,n);

	  printf("ZZ component:\n");
	  print_matrix(zz,n,n);

	  free(xx); free(xy); free(xz);
	  free(yy); free(yz); free(zz);
	  printf("\n</quadrupole>\n\n");
  }

  /*----------------------------------------------------------------
   * print one-electron Hamiltonian and kinetic energy integrals
   */
  if ( cap & HAVE_ONEEL ) {
	  int n=0;
	  double *h, *t;

	  phis_get_oneel(NULL, NULL, &n);
	  n = abs(n);

	  h = (double*)malloc_we(n*n*sizeof(double));
	  t = (double*)malloc_we(n*n*sizeof(double));
  
	  phis_get_oneel(h,t,&n);

	  printf("<!--------------------------------------------------------------"
		 "-------------->\n\n"
		 "<one-electron>\n\n");

	  printf("one-electron Hamiltonian:\n");
	  print_matrix(h,n,n);

	  printf("kinetic energy integrals:\n");
	  print_matrix(t,n,n);

	  free(h); free(t);
	  printf("\n</one-electron>\n\n");
  }

  /*----------------------------------------------------------------
   * print integrals 
   */
  if ( cap & HAVE_NEXTINT ) {
    printf("<!--------------------------------------------------------------"
	   "-------------->\n\n"
	   "<integrals>\n\n");

    if(print_integrals) {
      int p,q,r,s;
      double Vpqrs;
      
      printf("   p   q   r   s |  Vpqrs\n"
	     "-----------------+---------\n");

      phis_get_next_Vpqrs(&p,&q,&r,&s,&Vpqrs);
      do {
	printf(" %3d %3d %3d %3d | %7.4f\n",p,q,r,s,Vpqrs);
 	phis_get_next_Vpqrs(&p,&q,&r,&s,&Vpqrs);
      } while (s > 0);
    }
    else 
      printf("printing of two electron integrals skipped "
	     "(specify '-i' to see them)\n");

    printf("\n</integrals>\n\n");
  }

  return 0;
}


/*----------------------------------------------------------------*/


/* prints the matrices in a readable form */
void print_matrix(double * matrix,int nmb_rows,int nmb_cols)
{
  int i,j,k;
  for(k=0;k<(nmb_cols+9)/10;k++) {
    printf("     ");
    for(j=10*k;j<10*(k+1) && j<nmb_cols;j++) printf("|%10d  ",j+1);
    printf("\n");
    printf(" ----");
    for(j=10*k;j<10*(k+1) && j<nmb_cols;j++) printf("+------------");
    printf("\n");
    
    for(i=0;i<nmb_rows;i++) {
      printf(" %3d |", i+1);
      for(j=10*k;j<10*(k+1) && j<nmb_cols;j++) {
	printf("%12.8f ",matrix[nmb_rows*j+i]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

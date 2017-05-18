#include "phis.h"
#include "./guk.h"

void phis_guk_geometry(int *n,double *geo,double *Z_nuc,double *e)
{
  int n_words,i;
  double *tmp_geo;
  FILE *dump_fp;

  n_words=3*guk_maxat;

  /* check n, geo is expected to have size nx3 */
  if(*n<nAtoms) {
    *n=-nAtoms;
    return;
  }
  *n=nAtoms;

  /* allocate temporary space for geometry section */
  tmp_geo=(double*)malloc_we(n_words*sizeof(double));
  
  /* open dumpfile */
  if ((dump_fp = guk_open_dump_file(dumpfile_name)) == NULL) {
    perror( "Cannot open dumpfile" ); 
    exit( 1 );
  }
  
  /* read section type 5 (geometry information) */
  if(guk_read_data(dump_fp,(char*)tmp_geo,5,n_words)!=n_words) {
    printf("error reading zmatrix (section type 5) from dumpfile!\n");
    if(n_words) printf("%d words read.\n", n_words);
    perror("guk_read_data");
    exit(1);
  }

  /* copy geometry in array allocated by user */
  for(i=0;i<3*nAtoms;i++)
    geo[i]=tmp_geo[i];

  free(tmp_geo);


  /* copy nuclear charges from section type 1 */
  {
    Type_1 *section1;
    
    section1 = guk_get_section(dump_fp,1);
    for(i=0;i<nAtoms;i++)
      Z_nuc[i]=section1->z[i];

    free(section1);
  }


  /* nuclear energy is in section type 2 */
  guk_read_type2(e, 1, ENUC);
 
  fclose(dump_fp);
}

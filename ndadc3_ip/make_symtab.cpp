
#include <cmath>
extern "C" {
/*#include "phis.h"*/
#include "globals.h"
}
#include <iostream>
#include <string>
using namespace std;

//static char *cvsid="$Id: make_symtab.c,v 1.4 2003/12/01 14:54:55 joergb Exp $";

/* make_sysmtab:
 *
 * generates the symmetry tables for all following calculations
 * WARNING: In affinity mode the roles of particles and holes are swapped!
 *
 * $Log: make_symtab.c,v $
 * Revision 1.4  2003/12/01 14:54:55  joergb
 * parameters changed (make_symtab does not need the
 * whole inp structure => usable in other programs)
 *
 * Revision 1.3  2003/06/12 13:58:50  u85
 * now exists if input symmetry is empty
 *
 * Revision 1.2  2003/03/03 14:18:06  u85
 * Error message beautified
 *
 * Revision 1.1.1.1  2002/12/04 15:17:18  u85
 * Imported sources
 *
 */
extern "C" void make_symtab(int sym,int affinity_mode,int debug,FILE *data_fp,
			    struct scf_t *scf,struct symtab_t *symtab)
{
  int i;//,s,sym_index[9];
  char error[100];

  /* initialization */
  for(i=0;i<9;i++) {
    symtab->nOcc[i]=0;
    symtab->nVir[i]=0;
    symtab->occ[i]=0;
    symtab->vir[i]=0;
    //sym_index[i]=0;
  }

  /* check if symmetry specified in input is there */
  if(sym>scf->nSym) {
    sprintf(error,
            "Error: You want to calculate symmetry %d, but the number of\n"
            "symmetries in the molecule is just %d\n",
            sym,scf->nSym);
    throw string(error);
    //exit(42);
  }

  /* just copy number of symmetries tothis structure */
  symtab->nSym=scf->nSym;

  /* get the number of orbitals of each symmetry */
  for(i=1;i<=scf->nBas;i++) {
    if(fabs(scf->occ[i]-2.0)<SCF_THRESH)
      symtab->nOcc[scf->sym[i]]++;
    else if(scf->occ[i]<SCF_THRESH)
      symtab->nVir[scf->sym[i]]++;
    else {
         sprintf(error,
              "Error: Occupation of orbital %d is %5.2f , but your molecule\n"
              "has to be closed-shell\n",
              i,scf->occ[i]);
  
      throw string(error);
      //exit(42);
    }
  }

  /* allocate space for tabs */
  for(i=1;i<=symtab->nSym;i++) {
    symtab->occ[i]=new int [symtab->nOcc[i]+1];
    symtab->vir[i]=new int [symtab->nVir[i]+1];
    symtab->occ[i][symtab->nOcc[i]]=-1;
    symtab->vir[i][symtab->nVir[i]]=-1;
  }

  /* generate symmetry tabs */
  for(i=0;i<9;i++) {
    symtab->nOcc[i]=0;
    symtab->nVir[i]=0;
  }
  
  for(i=1;i<=scf->nBas;i++) {
    if(fabs(scf->occ[i]-2.0)<SCF_THRESH)
      symtab->occ[scf->sym[i]][symtab->nOcc[scf->sym[i]]++]=i;
    else if(scf->occ[i]<SCF_THRESH)
      symtab->vir[scf->sym[i]][symtab->nVir[scf->sym[i]]++]=i;
    else {
     sprintf(error,"occupation of orbital %d is %5.2f ,\n"
              "but your molecule has to be closed-shell\n",
              i,scf->occ[i]);

      throw string(error);
      //exit(42);
    }
  }
    
  /* swap virtual<->occupied in affinity mode */
  if(affinity_mode) {
    int dummy,*dummy_ptr;

    for(i=1;i<=8;i++) {
      dummy           = symtab->nOcc[i]; 
      symtab->nOcc[i] = symtab->nVir[i]; 
      symtab->nVir[i] = dummy;

      dummy_ptr       = symtab->occ[i]; 
      symtab->occ[i]  = symtab->vir[i]; 
      symtab->vir[i]  = dummy_ptr;
    }
  }  

  /* break if symmetry is empty */
  if(symtab->nOcc[sym]==0) {
    /* fprintf(stdout,"there is no orbital for a non-dyson calculation in " */
/*             "symmetry %d\n => nothing to be done\n", */
/*             sym); */
    /* exit(0);*/
    return;
  }

  /* debug output follows */

//   if(debug<2) return;
  
//   fprintf(data_fp,
// 	  "#------------------------------------------------------------"
// 	  "--------------\n\n"
// 	  "#begin<nmb_syms>\n\n"
// 	  "#irrep   occupied    virtual\n");

//   for(i=1;i<=symtab->nSym;i++)
//     fprintf(data_fp,
// 	    "  %2d       %4d       %4d\n",i,symtab->nOcc[i],symtab->nVir[i]);

//   fprintf(data_fp,
// 	  "\n#end<nmb_syms>\n\n");

//   fflush(data_fp);

//   if(debug<4) return;

//   fprintf(data_fp,
// 	  "#------------------------------------------------------------"
// 	  "--------------\n\n"
// 	  "#begin<symtab>\n\n"
// 	  "#irrep      label    orbital\n");

//   for(s=1;s<=symtab->nSym;s++) {
//     for(i=0;i<symtab->nOcc[s];i++)
//       fprintf(data_fp,
// 	      "  %2d       %4d*      %4d\n",
// 	      s,i,symtab->occ[s][i]);
//     for(i=0;i<symtab->nVir[s];i++)
//       fprintf(data_fp,
// 	      "  %2d       %4d       %4d\n",
// 	      s,i,symtab->vir[s][i]);
//   }

//   fprintf(data_fp,
// 	  "\n#end<symtab>\n\n");

//   fflush(data_fp);
}

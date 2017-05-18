#include <stdlib.h>

#include "my_globals.h"
#include "filetools.h"
#include "my_adc_macros.h"
#include "read_scf_data.h"
#include <math.h>
#include <string.h>

#define KL (*k*(*k-1)/2+*l-1)

static char *cvsid="$Id: calc_psi0_rho_psi0.c,v 1.2 2004/05/09 17:10:00 imke Exp $";

/* cap_calc_psi0_rho_psi0:
 *
 * calculates ground state expectation value of rho
 *
 * $Log: calc_psi0_rho_psi0.c,v $
 * Revision 1.2 2004/05/09 17:10 imke
 *
 * Revision 1.1  2004/01/29 16:57:26  joergb
 * Initial import
 *
 */
void cap_calc_psi0_rho_psi0(struct inp_t *inp,struct scf_t *scf,
			struct symtab_t *symtab,double *rho)
{
  int *k,*l,*sym,i;
  double epsi_k,epsi_l;
  double *sigma,*q_kl;
  char *sigma_file=(char*)inp->sigma_file;
  char *q_kl_file=(char*)inp->q_kl_file;

  printf("calculating < psi_0 | rho | psi_0 >\n"
	 "(using sigma and q_kl files)\n\n");


  /* rho = delta_kl n_k = < phi_0 | c_k^+ c_l | phi_0 > */
  /* stimmt das? */
  for(i=1;i<=scf->nBas;i++) 
    if(scf->occ[i]>2.-SCF_THRESH) 
      rho[i*(i-1)/2+i-1]=1.;
    

  for(sym=inp->sigma_file_type+1;*sym!=-1;sym++) {

    printf("Initialize file of symmtery %i\n", *sym);

    /* initialize files */
    inp->sigma_file=init_fileinfo(sigma_file,inp->sigma_file_type,
				  -1,*sym,RD_MODE);
    inp->q_kl_file =init_fileinfo(q_kl_file, inp->sigma_file_type,
				  -1,*sym,RD_MODE);


    /* allocate matrices */
    make_matrix(inp->sigma_file->dim,inp->sigma_file->dim,'p','n',&sigma);
    make_matrix(inp->sigma_file->dim,inp->sigma_file->dim,'p','n',&q_kl);



    /* read matrices */
    FT_READ_MATRIX(inp->sigma_file,sigma);
    printf(" finished read sigma_file of sym %i\n", *sym);
    FT_READ_MATRIX(inp->q_kl_file ,q_kl);
    printf(" finished reading qkl-file of sym %i\n", *sym);


    /* rho += q_kl */
    FOR_OCC_IN_SYM(k,*sym) {
      FOR_OCC_IN_SYM(l,*sym) {
	if(l>k) break;

	rho[KL] += q_kl[XO(k)*(XO(k)+1)/2+XO(l)];
      }
    }

    /* new rho += q_kl fuer alle! */
    FOR_VIR_IN_SYM(k,*sym){
      FOR_OCC_IN_SYM(l,*sym){
	if(l>k) break;

	rho[KL] +=q_kl[(XV(k))*(XV(k)+1)/2+XO(l)];
      }
    }

    FOR_VIR_IN_SYM(k,*sym){
      FOR_VIR_IN_SYM(l,*sym){
	if(l>k) break;

	rho[KL] +=q_kl[(XV(k))*(XV(k)+1)/2+XV(l)];
      }
    }

    /* Ende der Aenderung */


    /* rho += -{\bar n}_k n_l/(epsi_k-epsi_l) * sigma_lk  */
    FOR_VIR_IN_SYM_WE(k,*sym,epsi_k) { 
      FOR_OCC_IN_SYM_WE(l,*sym,epsi_l) {
	rho[KL] += -sigma[(XV(k))*(XV(k)+1)/2+XO(l)]/(epsi_k-epsi_l);
      }
    }


    /* free memory */
    free(sigma);
    free(q_kl);


    /* close files */
    /*    FT_CLOSE(inp->sigma_file);
	  FT_CLOSE(inp->q_kl_file);*/

    printf("\n");
  }
}


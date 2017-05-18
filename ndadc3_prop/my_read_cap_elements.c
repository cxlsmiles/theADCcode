#include <stdio.h>
#include <stdlib.h>

#include "phis.h"
#include "my_globals.h"

/* read_cap_elements: */
/*                                                                          */
/* reads the ASCII output of the secap program (outform=1)                  */
/* theses files (one for each non-empty irrep) must be present named        */
/* named Wse.number_of_irrep                                                */
/*                                                                          */
/* first: check if number of virtual is consistent to symtab                */
/* reads first diagonal elements                                            */
/* then off-diagonal elements of vitual/virtual combinations only           */
/* storage is nevertheless in a upper-triangonal by-colum matrix including  */
/* all 0's for the occupied orbitals because the d_matrix is supposed to be */
/* calculated for any given one-particle operator read from this format     */
/* furthermore, all 0's between different irreps are also stored ...        */
            

/* please notice: only molcas is supported by secap so far, this may change */
/* if secap is based of phis one day, output stayes hoepfully the same :-)  */

/* notice: virtuals are in energetical ordering, NOT in symmetry ordering!! */
/* global indices must thus be used throughout this module                  */


void read_cap_elements(struct inp_t *inp, struct symtab_t *symtab, double *cap_elements)
{
  FILE *cap_file_ptr;
  int i=0, j=0, k=0; 

  int index,index_i,index_j;

  int cap_virt;

  int num_virt=0, num_occ=0;          /* 0's are neccessary here! */
  int num_orb;

  char buch[11];                          /* helper for fscanf */

  char buf[9][10];                   /* buffer for reading the filenames */

  /*----------------------------------------------------------------------------------*/
  /* calculate number of virtual and occupied orbitals from symtab                    */
  /* seems not to be used, delete....                                                 */
  /* (was used for allocating cap_elements, which is now done in calling routine ...) */
  /*----------------------------------------------------------------------------------*/
  for(i=1; i<= (symtab->nSym); i++){
    num_virt=num_virt + symtab->nVir[i];
    num_occ =num_occ  + symtab->nOcc[i];
  }
  num_orb= num_occ + num_virt;

  /*--------------------------------------------------------------------------*/
  /* open cap_file with name Wse.numberofirrep if irrep not empty .....       */
  /* and check if file was opened ....                                        */

  for(k=1; k<=symtab->nSym; k++){
    printf(" number of virtuals: %i \n", symtab->nVir[k]); 
    if(symtab->nVir[k]!=0){
      sprintf(buf[k],"%s.%i",inp->cap_file, k);
      cap_file_ptr=fopen(buf[k], "r");
      printf("open File %s ...\n", buf[k]);

      if(!cap_file_ptr){
	printf("failed open %s, please ensure that %s file is present in work directory\n", buf[k], buf[k]);
        printf("program will stop.\n");
	exit(23);
      }

      /*---------------------------------------------------------------------------*/
      /* read cap_file and compare number of virtuals to symtab number of virtuals */
      /* overread char virtuals and enmpty line in Wse.1                           */
      /* notice: Cap_file with name Wse.1 expected!                                */
      /*---------------------------------------------------------------------------*/
      fscanf(cap_file_ptr, "%i", &cap_virt);
      fscanf(cap_file_ptr, "%s\n", buch);

      if(cap_virt!=symtab->nVir[k]){
	printf("Number of virtuals in cap_file ( %i ) does not match number of virtuals \n read from symtab ( %i ) \n program will stop.\n", cap_virt, symtab->nVir[k]);
	exit(24);
      }

      /*---------------------------------------------------------------------------*/
      /* allocate matrix cap_elements from symtab data:                            */
      /* initialize the matrix cap_elements with 0.0's                             */
      /* (FIXME: might be that allocation has to be done in the main program       */
      /* if so, just transfer these lines there.)                                  */
      /* ANYWAY: cap_elements has to be global!!                                   */
      /* already done by make_matrix in calling routine                            */
      /*---------------------------------------------------------------------------*/

      /*cap_elements=(double*)malloc(sizeof(double)*(num_orb*(num_orb+1)/2));
      for(i=0; i < (num_orb*(num_orb+1))/2 ; i++){
	cap_elements[i]=0.0;
	}*/

      /*--------------------------------------------------------------------------*/
      /* read virtual diagonal elements from cap_file and store them on position :*/
      /* global_index * (global_index +1) /2 -1                                   */
      /* notice: c-numbering starts at 0 (cap_elements)!!                         */
      /* (FIXME: check if i==index_i-1)                                           */ 
      /* do not read empty line in Wse.x file                                     */
      /*--------------------------------------------------------------------------*/

      for(i=0; i< symtab->nVir[k]; i++){   /* Speicherung in symtab geht mit 0 los */
    	index=((symtab->vir[k][i])*((symtab->vir[k][i]) +1))/2 -1;
	fscanf(cap_file_ptr, "%i %lf", &index_i, &cap_elements[index]);
	if(inp->debug > 0){
	  printf(" irrep: %i symtab_index: %i read_index: %i matrix_index: %i element: %16.14f\n", k, symtab->vir[k][i], index_i, index, cap_elements[index]);
	}
      }

      fscanf(cap_file_ptr, "\n");             /* do not read empty line */


      /*--------------------------------------------------------------*/
      /* read off_diagonal matrix elements and store them on position:*/
      /* ((global_index_i*(global_index_i -1))/2 -1 + global_index_j  */
      /*  |||||||||||||||||||||||||||||||||||||||||                   */
      /*   diag element of previous col.                              */
      /*                                              ||||||||||||    */
      /*                                       posit. in present col. */
      /* (FIXME: check if i==index_i -1 and j==index_j -1 )           */
      /* note: j<i is asserted by this loop and a must!!              */
      /*--------------------------------------------------------------*/

      for(i=1; i< symtab->nVir[k]; i++){     /* Speicherung beginnt bei 0, d.h. 1 heisst */
	for(j=0; j<i; j++){              /* virtuelles orb Nr. 2, ebenso j=0, virt Nr.1 */
	  index=((symtab->vir[k][i])*((symtab->vir[k][i]) -1))/2 -1 + (symtab->vir[k][j]);
	  fscanf(cap_file_ptr, "%i %i %lf", &index_i, &index_j, &cap_elements[index]);
	  if(inp->debug > 2){
	    printf(" irrep: %i symtab_index: %i %i read_indeces: %i %i matrix_index: %i element: %16.14f\n", k, symtab->vir[k][i], symtab->vir[k][j], index_i, index_j, index, cap_elements[index]);
	  }
	}
      }

  /*-------------------------------------------------------------------------*/
  /* That's it, we are done  for one irrep                                   */
  /* close file (a Must because only one filepointer is used!), some output! */
  /*-------------------------------------------------------------------------*/
      fclose(cap_file_ptr);
      printf("finished reading cap_elements from %s: %i virtuals have been identified\n", buf[k], cap_virt);
    }
  }

  /* debug ausgabe */
  if(inp->debug >3){
    for(i = 0; i<(num_orb*(num_orb+1))/2; i++){
      printf("Globalindex: %i Element: %16.14f\n", i, cap_elements[i]);
    }
  }

  /* all elements should be read by now ... */
  printf("finished reading (MO/CAP/MO) for all irreps with virtual orbitals\n");
  printf("leaving routine read_cap_elements, bye!\n");

}



/*---------------------------------------------------------------------------------------------*/
/* void read_all_cap_elements(struct inp_t *inp, struct symtab_t *symtab, double *cap_elements)*/
/* like read_cap_elements, but reads (MO/CAP/MO) for ALL orbitals occupied/occupied,           */
/* unoccupied/occupied  and unoccupied/unoccupied                                              */


/* read_all_cap_elements: */
/*                                                                          */
/* reads the ASCII output of the secap program (outform=1) (allorbs=1)      */
/* theses files (one for each non-empty irrep) must be present named        */
/* named Wse.number_of_irrep                                                */
/*                                                                          */
/* first: check if number of orbitals is consistent to symtab               */
/* reads first diagonal elements                                            */
/* then off-diagonal elements of MO/MO combinations                         */
/* storage in a upper-triangonal by-colum matrix                            */
            

/* please notice: only molcas is supported by secap so far, this may change */
/* if secap is based of phis one day, output stayes hoepfully the same :-)  */

/* notice: MOs are in energetical ordering, NOT in symmetry ordering!!      */
/* global indices must thus be used throughout this module                  */



void read_all_cap_elements(struct inp_t *inp, struct symtab_t *symtab, double *cap_elements)
{
  FILE *cap_file_ptr;
  int i=0, j=0, k=0; 
  int i_help, j_help;

  int index,index_i,index_j;

  int cap_mo;

  int num_virt=0, num_occ=0, numorb=0;          /* 0's are neccessary here! */
  int num_orb_all;

  char buch[11];                          /* helper for fscanf */

  char buf[9][10];                   /* buffer for reading the filenames */

  /*----------------------------------------------------------------------------------*/
  /* calculate overall number of virtual and occupied orbitals from symtab            */
  /*----------------------------------------------------------------------------------*/
  for(i=1; i<= (symtab->nSym); i++){
    num_virt=num_virt + symtab->nVir[i];
    num_occ =num_occ  + symtab->nOcc[i];
  }
  num_orb_all= num_occ + num_virt;

  /*--------------------------------------------------------------------------*/
  /* calc number of MOs for current irrep                                     */
  /* open cap_file with name Wse.numberofirrep if irrep not empty .....       */
  /* and check if file was opened ....                                        */

  for(k=1; k<=symtab->nSym; k++){
    numorb=symtab->nVir[k]+symtab->nOcc[k];
    printf(" number of MOs: %i   (number of vituals: %i number of occupied: %i)\n", numorb, symtab->nVir[k], symtab->nOcc[k]); 
    if(numorb!=0){
      sprintf(buf[k],"%s.%i",inp->cap_file, k);
      cap_file_ptr=fopen(buf[k], "r");
      printf("open File %s ...\n", buf[k]);

      if(!cap_file_ptr){
	printf("failed open %s, please ensure that %s file is present in work directory\n", buf[k], buf[k]);
        printf("program will stop.\n");
	exit(23);
      }

      /*---------------------------------------------------------------------------*/
      /* read cap_file and compare number of virtuals to symtab number of virtuals */
      /* overread char virtuals and enmpty line in Wse.1                           */
      /* notice: Cap_file with name Wse.1 expected!                                */
      /*---------------------------------------------------------------------------*/
      fscanf(cap_file_ptr, "%i", &cap_mo);
      fscanf(cap_file_ptr, "%s\n", buch);

      if(cap_mo!=numorb){
	printf("Number of MOs in cap_file ( %i ) does not match number of MOs \n read from symtab ( %i ) \n program will stop.\n", cap_mo, numorb);
	exit(24);
      }

      /*--------------------------------------------------------------------------*/
      /* read  diagonal elements from cap_file and store them on position :       */
      /* global_index * (global_index +1) /2 -1                                   */
      /* notice: first symtab->nOcc[k] orbitals belong to occupied, rest to virt! */
      /* notice: c-numbering starts at 0 (cap_elements)!!                         */
      /* (FIXME: check if i==index_i-1)                                           */ 
      /* do not read empty line in Wse.x file                                     */
      /*--------------------------------------------------------------------------*/

      for(i=0; i< numorb; i++){   /* Speicherung in symtab geht mit 0 los */
	if(i < symtab->nOcc[k]){
	  index=((symtab->occ[k][i])*((symtab->occ[k][i]) +1))/2 -1;
	  fscanf(cap_file_ptr, "%i %lf", &index_i, &cap_elements[index]);
	  if(inp->debug > 0){
	    printf(" irrep: %i occ. symtab_index: %i read_index: %i matrix_index: %i element: %16.14f\n", k, symtab->occ[k][i], index_i, index, cap_elements[index]);
	  }
	}
	else{
	  i_help=i-(symtab->nOcc[k]);
	  index=((symtab->vir[k][i_help])*((symtab->vir[k][i_help]) +1))/2 -1;
	  fscanf(cap_file_ptr, "%i %lf", &index_i, &cap_elements[index]);
	  if(inp->debug > 0){
	    printf(" irrep: %i virt. symtab_index: %i read_index: %i matrix_index: %i element: %16.14f\n", k, symtab->vir[k][i_help], index_i, index, cap_elements[index]);
	  }
	}
      }

      fscanf(cap_file_ptr, "\n");             /* do not read empty line */


      /*--------------------------------------------------------------------*/
      /* read off_diagonal matrix elements and store them on position:      */
      /* ((global_index_i*(global_index_i -1))/2 -1 + global_index_j        */
      /*  |||||||||||||||||||||||||||||||||||||||||                         */
      /*   diag element of previous col.                                    */
      /*                                              ||||||||||||          */
      /*                                       posit. in present col.       */
      /* (FIXME: check if i,index_i -1 consitent and j,index_j -1 cons.     */
      /* note: j<i is asserted by this loop and a must!!                    */
      /* also a must: check if global index has to be read from virt or occ */
      /* -> i_help, j_help  (= i-number of occupied )                       */
      /*--------------------------------------------------------------------*/

      for(i=1; i< numorb; i++){     /* Speicherung beginnt bei 0, d.h. 1 heisst */
	for(j=0; j<i; j++){              /*  orb Nr. 2, ebenso j=0, MO Nr.1     */
	  if(i< symtab->nOcc[k]){
	    index=((symtab->occ[k][i])*((symtab->occ[k][i]) -1))/2 -1 + (symtab->occ[k][j]);
	    fscanf(cap_file_ptr, "%i %i %lf", &index_i, &index_j, &cap_elements[index]);
	    if(inp->debug > 2){
	      printf(" irrep: %i symtab_index (occ./occ.): %i %i read_indeces: %i %i matrix_index: %i element: %16.14f\n", k, symtab->occ[k][i], symtab->occ[k][j], index_i, index_j, index, cap_elements[index]);
	    }
	  }
	  else{
	    i_help=i-(symtab->nOcc[k]);
	    if(j< symtab->nOcc[k]){
	      index=((symtab->vir[k][i_help])*((symtab->vir[k][i_help]) -1))/2 -1 + (symtab->occ[k][j]);
	      fscanf(cap_file_ptr, "%i %i %lf", &index_i, &index_j, &cap_elements[index]);
	      if(inp->debug > 2){
	      printf(" irrep: %i symtab_index (vir./occ.): %i %i read_indeces: %i %i matrix_index: %i element: %16.14f\n", k, symtab->vir[k][i_help], symtab->occ[k][j], index_i, index_j, index, cap_elements[index]);
	      }
	    }
	    else{
	      j_help=j-(symtab->nOcc[k]);
	      index=((symtab->vir[k][i_help])*((symtab->vir[k][i_help]) -1))/2 -1 + (symtab->vir[k][j_help]);
	      fscanf(cap_file_ptr, "%i %i %lf", &index_i, &index_j, &cap_elements[index]);
	      if(inp->debug > 2){
		printf(" irrep: %i symtab_index (vir./vir.): %i %i read_indeces: %i %i matrix_index: %i element: %16.14f\n", k, symtab->vir[k][i_help], symtab->vir[k][j_help], index_i, index_j, index, cap_elements[index]);
	      }
	    }
	  }
	}
      }


  /*-------------------------------------------------------------------------*/
  /* That's it, we are done  for one irrep                                   */
  /* close file (a Must because only one filepointer is used!), some output! */
  /*-------------------------------------------------------------------------*/
      fclose(cap_file_ptr);
      printf("finished reading cap_elements from %s: %i virtuals have been identified\n", buf[k], cap_mo);
    }
  }

  /* debug ausgabe */
  if(inp->debug >3){
    for(i = 0; i<(num_orb_all*(num_orb_all+1))/2; i++){
      printf("Globalindex: %i Element: %16.14f\n", i, cap_elements[i]);
    }
  }

  /* all elements should be read by now ... */
  printf("finished reading (MO/CAP/MO) for all irreps with active MOs \n");
  printf("leaving routine read_all_cap_elements, bye!\n");

}


#include <string.h>

#include "./guk.h"

Type_1*
guk_init_Type_1(Type_1* sec)
{
	/* if called with NULL pointer allocate structure */
	if(!sec) {
		sec=(Type_1*)malloc(sizeof(Type_1));
		sec->data=NULL;
	}

	/* set pointers and determine size of data area */
	sec->ex         = (double*)(sec->data);
	sec->cc         = (double*)(sec->ex          + guk_mxprim);       
	sec->z          = (double*)(sec->cc          + guk_maxaqm*guk_mxprim);
	sec->job_title  =   (char*)(sec->z           + guk_maxat);        
	sec->atom_labels=   (char*)(sec->job_title   + 80);           
	sec->kstart     =    (int*)(sec->atom_labels + guk_maxat*8);      
	sec->katom      =    (int*)(sec->kstart      + guk_mxshel);       
	sec->ktype      =    (int*)(sec->katom       + guk_mxshel);       
	sec->kng        =    (int*)(sec->ktype       + guk_mxshel);       
	sec->kloc       =    (int*)(sec->kng         + guk_mxshel);       
	sec->kmin       =    (int*)(sec->kloc        + guk_mxshel);       
	sec->kmax       =    (int*)(sec->kmin        + guk_mxshel);       
	sec->size       =   (char*)(sec->kmax        + guk_mxshel + 3) - sec->data;

	if(sec->data) {
		sec->n_shells = sec->kmax[guk_mxshel];             
		sec->n_atoms  = sec->kmax[guk_mxshel+1];             
		sec->n_basis  = sec->kmax[guk_mxshel+2];             
	}
	else
		/* allocate space for data (stored in words => "+8") */
		sec->data=(char*)malloc(sec->size+8);

	return(sec);
}


Type_3*
guk_init_Type_3(Type_3* sec)
{
	/* if called with NULL pointer allocate structure */
	if(!sec) {
		sec=(Type_3*)malloc(sizeof(Type_3));
		sec->data=NULL;
	}

	/* set pointers and determine size of data area */
	sec->user      =   (char*)(sec->data);
	sec->date      =   (char*)(sec->user      + 8);
	sec->time      =   (char*)(sec->date      + 8);
	sec->prog      =   (char*)(sec->time      + 8);
	sec->scftype   =   (char*)(sec->prog      + 8);
	sec->acct      =   (char*)(sec->scftype   + 8);
	sec->more      =   (char*)(sec->acct      + 8);
	sec->job_title =   (char*)(sec->more      + 104);
	sec->orb_eng   = (double*)(sec->job_title + 80);
	sec->occ_num   = (double*)(sec->orb_eng   + guk_maxorb);
        /* 1 double + 3 ints */
	sec->unused    =    (int*)(sec->occ_num   + guk_maxorb + 1 ) + 3;
	sec->ilifc     =    (int*)(sec->unused    + 3);
	sec->ntran     =    (int*)(sec->ilifc     + guk_maxorb);
	sec->itran     =    (int*)(sec->ntran     + guk_maxorb);
	sec->ctran     = (double*)(sec->itran     + guk_maxorb*3);
	/* 1 int */
	sec->size      =  ((char*)(sec->ctran     + guk_maxorb*3) - sec->data) + sizeof(int);

	if(sec->data) {
		sec->ehf_tot = sec->occ_num[guk_maxorb];
		sec->n_gtos  = sec->unused[-3];
		sec->nlcbf   = sec->unused[-2];
		sec->n_vec   = sec->unused[-1];
		sec->iftran  = *(int*)(sec->ctran+guk_maxorb*3);
	}
	else
		/* allocate space for data (stored in words => "+8") */
		sec->data = (char*)malloc(sec->size+8);

	return(sec);
}


Type_51*
guk_init_Type_51(Type_51* sec)
{
	/* if called with NULL pointer allocate structure */
	if(!sec) {
		sec=(Type_51*)malloc(sizeof(Type_51));
		sec->data=NULL;
	}

	/* set pointers and determine size of data area */
	/* 1 int */
	sec->mult_table   =  (int*)(sec->data)        + 1;
	sec->symlabels_ao =  (int*)(sec->mult_table   + 64);
	sec->symlabels_mo =  (int*)(sec->symlabels_ao + guk_maxorb);
	sec->size         = (char*)(sec->symlabels_mo + guk_maxorb) - sec->data;

	if(sec->data) {
		sec->n_irrep = sec->mult_table[-1];
	}
	else
		/* allocate space for data (stored in words => "+8") */
		sec->data = (char*)malloc(sec->size+8);

	return(sec);
}


Type_1005*
guk_init_Type_1005(Type_1005* sec)
{
	/* if called with NULL pointer allocate structure */
	if(!sec) {
		sec=(Type_1005*)malloc(sizeof(Type_1005));
		sec->data=NULL;
	}

	/* set pointers and determine size of data area */
	/* 1 int */
	sec->list_of_act =  (int*)(sec->data)       + 1;
	/* 2 ints */
	sec->list_of_fro =  (int*)(sec->list_of_act + guk_maxorb + 2);
	/* 1 int */
	sec->size        = (char*)(sec->list_of_fro + guk_maxorb + 2) - sec->data;

	if(sec->data) {
		sec->nAct           = sec->list_of_act[-1];
		sec->act_print_flag = sec->list_of_act[guk_maxorb];
		sec->nFro           = sec->list_of_fro[-1];
		sec->fro_print_flag = sec->list_of_fro[guk_maxorb];
	}
	else
		/* allocate space for data (stored in words => "+8") */
		sec->data = (char*)malloc(sec->size+8);

	return(sec);
}


char*
guk_init_section(char* section, int type)
{
	switch (type) {
	case 1 :
		section = (char*)guk_init_Type_1((Type_1*) section);
		break;
	case 3 :
		section = (char*)guk_init_Type_3((Type_3*) section);
		break;
	case 51 :
		section = (char*)guk_init_Type_51((Type_51*) section);
		break;
	case 1005 :
		section = (char*)guk_init_Type_1005((Type_1005*) section);
		break;
	default :

		printf("Sections of type %i are not supported.\n", type);
		return NULL;
		break; 
	}

	return(section);
}


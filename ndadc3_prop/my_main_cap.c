static char *cvsid="$Id: main_cap.c,v 1.0 2004/ /  time imke Exp $";



/* main_cap.c                                             */
/* main source of program generating the CAP-matrix       */
/* for non-dyson ADC(3) using the ISR                     */
/*                                                        */
/* required: integrals, scf, transformation               */
/* of molecule, only molcas5 backend supported until      */
/* external cap-mo-representation is based on phis        */
/* or got a backend for another program package           */
/*                                                        */
/* required: file Wse.1, ASCII file containing the        */
/* (MO/CAP/MO) integrals provided that both MO are        */
/* virtuals.                                              */
/* check: use secap5 for molcas 5, don't forget to create */
/* ana.log file from seward(?) output                     */
/*                                                        */
/*                                                        */
/* 1) get input using el_parso input parser, defaults     */ 
/*    are set.                                            */
/* 2) read the cap-representation from Wse.1              */
/* 3) build the ISR of the cap                            */
/* 4) write the resulting matrix as upper triangel by col */
/*    to file specified in input.                         */
/*                                                        */
/* $Log: main_cap.c, v $                                  */
/* Revision 1.0  date time imke                           */
/* Initial import                                         */
/*                                                        */
/*--------------------------------------------------------*/


#include "my_globals.h"
#include <stdlib.h>                                                        
#include "jbutils.h"
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{

  struct inp_t inp;
  struct scf_t scf;
  struct symtab_t symtab;
  char *defaults;

  int parse_stdin=1,c;
  char hostname[20];

  /*-------------------------------------------------------*/
  /*  default input value for run based on molcas backend  */
  /*  reminder: use condig for generation of sigma and qkl!*/
  /*-------------------------------------------------------*/

  char default_input[]=
       "backend     =\"molcas\";"
       "debug       =0;"
       "cap_file    =\"Wse\";"
       "cap_file_type =0;"
       "dout_file_type=1,0;"   /* new!*/
       ;

  /*------------------------------------------------------*/
  /* some initial output, reminder of molcas compatibility*/
  /*------------------------------------------------------*/

  gethostname(hostname,20);
  printf("\n            ************************************"
         "\n              cap generation module for nd_adc3"
         "\n            ************************************"
         "\n                                              \n"
         "\n remember: only molcas backend supported at present"
         "\n please condig to create sigma and qkl and reorder"
         "\n the irreps according to nd_adc3 ordering!");

  defaults=default_input;      
  parse_stdin=1;

  /*-----------------------------------------------------*/
  /* command line parsing by getopt                      */
  /*-----------------------------------------------------*/
  do{
    c = getopt(argc,argv,"h");
    switch(c){
    case -1: break;
    case 'h':
      defaults="help_all;";
      parse_stdin=-1;
      c = -1;
      break;
    case '?':
      fflush(stdout);
      fprintf(stderr,"('-h' for help)\n");
      exit(42);
    default:
      fflush(stdout);
      fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
      exit(42);
    }
  }while(c!=-1);


  /*-----------------------------------------------------------*/
  /* parse defaults and stdin (if required)                    */
  /*-----------------------------------------------------------*/

  get_input_cap(defaults,parse_stdin,&inp);
  printf("finished parsing input\n");

  /*-----------------------------------------------------------*/
  /* call matrix_builder: all real work is done here to keep   */
  /* this file short!                                          */
  /*-----------------------------------------------------------*/

  printf("starting building matrices\n");
  build_d_matrices(&inp,&scf,&symtab);

  return 0;
}





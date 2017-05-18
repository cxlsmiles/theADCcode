#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

void *malloc_we(size_t size)
{
  void *ret;

  if((ret=malloc(size))!=NULL) 
    return(ret);

  /* here, if error occured */

  fprintf(stderr,"problem allocating memory!\n");
  perror("malloc");
  exit(-1);
}

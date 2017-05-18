CPP      = /cvos/shared/apps/gcc/4.3.2/bin/cpp
CPPFLAGS = 

CC     = gcc
CFLAGS = -I../include -g -O2 

FC     = pgf90
FFLAGS = -tp k8-64 -fastsse -Mipa=fast 

OBJS= ../lib/libphis.a(init.o Vpqrs.o version.o util.o \
blas/dgemm.o blas/lsame.o  blas/xerbla.o)
cleanfiles     = *~ *.o 
realcleanfiles = $(cleanfiles) Makefile
distcleanfiles = $(realcleanfiles)

.PHONY: test

all: $(OBJS)

# Make sure to rebuild all relevant files if the configuration
# changes (automake would take care of this, too).
Vpqrs.o: ../include/phis.h ../include/config.h
init.o: ../include/phis.h ../include/config.h
version.o: ../include/phis.h ../include/config.h
util.o:

test:

include ../include/cleanup.mk

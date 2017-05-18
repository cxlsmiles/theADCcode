CPP      = /usr/bin/cpp
CPPFLAGS = 

CC     = icc
CFLAGS = -I../include -axT -xT -O3 -ip 

FC     = ifort
FFLAGS = -axT -xT -O3 -ip -pad 

OBJS= ../lib/libphis_guk.a(init.o Vpqrs.o version.o util.o \
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

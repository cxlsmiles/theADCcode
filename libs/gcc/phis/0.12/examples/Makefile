CC     = gcc
CFLAGS = -I../include -g -O2 

FC     = gfortran
FFLAGS = -fno-second-underscore -O2 -funroll-loops -fomit-frame-pointer  -mtune=opteron -mfpmath=sse,387 -m64 -march=opteron -mcmodel=medium

LDFLAGS = 
LIBS    = ../lib/libphis.a  /lib/libmolcas.a /libma.a
FLIBS   =  -L/cvos/shared/apps/gcc/4.3.2/lib/gcc/x86_64-unknown-linux-gnu/4.3.2 -L/cvos/shared/apps/gcc/4.3.2/lib/gcc/x86_64-unknown-linux-gnu/4.3.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/cvos/shared/apps/gcc/4.3.2/lib/gcc/x86_64-unknown-linux-gnu/4.3.2/../../.. -lgfortranbegin -lgfortran -lm

cleanfiles     = *.o *~ $(all)
realcleanfiles = Makefile sample fsample dump_all fdump_all int_check fint_check
distcleanfiles =

.PHONY: test

all: sample fsample dump_all fdump_all int_check fint_check

dump_all: dump_all.o ../lib/libphis.a
	$(CC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

int_check: int_check.o ../lib/libphis.a
	$(CC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

fint_check: fint_check.o ../lib/libphis.a
	$(FC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

fdump_all: fdump_all.o ../lib/libphis.a
	$(FC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

sample: sample.o ../lib/libphis.a
	$(CC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

fsample: fsample.o ../lib/libphis.a
	$(FC) $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

test:

include ../include/cleanup.mk

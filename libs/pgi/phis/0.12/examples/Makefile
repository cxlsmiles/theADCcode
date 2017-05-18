CC     = gcc
CFLAGS = -I../include -g -O2 

FC     = pgf90
FFLAGS = -tp k8-64 -fastsse -Mipa=fast 

LDFLAGS = 
LIBS    = ../lib/libphis.a  /lib/libmolcas.a /libma.a
FLIBS   =  -L/cvos/shared/TC/pgi/10.9/linux86-64/10.9/lib -L/usr/lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lnspgc -lpgc -lrt -lpthread -lm

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

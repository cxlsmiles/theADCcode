CC     = icc
CFLAGS = -I../include -axT -xT -O3 -ip 

FC     = ifort
FFLAGS = -axT -xT -O3 -ip -pad 

LDFLAGS = 
LIBS    = ../lib/libphis.a  /lib/libmolcas.a /libma.a
FLIBS   =  -L/cvos/shared/apps/intel/Compiler/11.1.046/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../.. -L/lib64 -L/lib -L/usr/lib64 -L/usr/lib -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lpthread -lirc_s -ldl

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

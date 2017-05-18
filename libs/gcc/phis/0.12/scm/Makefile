CC     = gcc
CFLAGS = -I../include -g -Wall

FC     = ifc
FFLAGS = -g

LDFLAGS =  -L/usr/local/molcas5/lib -L/usr/local/molcas5/g/lib/LINUX
LIBS    = ../lib/libphis.a -lmolcas -lma /usr/lib/gcc-lib/i486-linux/3.3.2/libg2c.a -lefence
# FLIBS   =  -L/usr/local/intel/compiler70/ia32/lib -L/usr/lib -lintrins -lCEPCF90 -lF90 -limf -lm -lirc -lcxa -lunwind

phis-scm.so: phis-scm.o ../lib/libphis.a
	$(CC) -shared $(LDFLAGS) $< $(LIBS) $(FLIBS) -o $@

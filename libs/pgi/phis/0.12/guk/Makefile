CC     = gcc
CFLAGS   = -I../include -g -O2

CPP      = /cvos/shared/apps/gcc/4.3.2/bin/cpp
CPPFLAGS = 

LDFLAGS  = 
LIBS     = ../lib/libphis.a  /lib/libmolcas.a /libma.a

RM       = /bin/rm
PERL     = /usr/local/bin/perl

# dependencies
OBJS= ../lib/libphis_guk.a(guk_init.o guk_aux.o sort_gukscf.o \
typeinits.o guk_read_type2.o \
phis_guk_info.o phis_guk_epsi.o phis_guk_sym.o phis_guk_occ.o \
phis_guk_nextint.o phis_guk_scfvec.o phis_guk_overlap.o \
phis_guk_ao.o phis_guk_loa.o phis_guk_geometry.o phis_guk_dip.o)

cleanfiles     = *.o *~
realcleanfiles = Makefile
distcleanfiles = 
.PHONY: test

# actions
all: $(OBJS)

test: 
	cd test/; \
	$(PERL) test.pl

include ../include/cleanup.mk

RM = /bin/rm

SUBDIRS = phis guk

cleanfiles     = *~ *.o lib/libphis_guk.a lib/libphis_guk.so
realcleanfiles = config.log config.status Makefile include/config.h autom4te.cache ./include/molcas.fh
distcleanfiles = configure ./include/molcas.fh

all:
	test -e lib || mkdir lib;
	@for subdir in $(SUBDIRS); do \
	  ( cd $$subdir && $(MAKE) $@ ) || exit 1; \
	done;
	ifort -shared  `ls */.libs/*.o` -o ./lib/libphis_guk.so
	rm `ls */*.lo`

test:
	@for subdir in $(SUBDIRS); do \
	  ( cd $$subdir && $(MAKE) $@ ) || exit 1; \
	done;

clean:
	@for subdir in $(SUBDIRS); do \
	  ( cd $$subdir && $(MAKE) $@ ) || exit 1; \
	done;
	-rm -rf $(cleanfiles)

# after this you will need to reconfigure and recompile
realclean: 
	@for subdir in $(SUBDIRS); do \
	  ( cd $$subdir && $(MAKE) $@ ) || exit 1; \
	done;
	-rm -rf $(cleanfiles) $(realcleanfiles)
	-rm -rf lib

# removes everything that is not in the CVS
distclean:
	@for subdir in $(SUBDIRS); do \
	  ( cd $$subdir && $(MAKE) $@ ) || exit 1; \
	done;
	-rm -rf $(cleanfiles) $(realcleanfiles) $(distcleanfiles)
	-rm -rf lib

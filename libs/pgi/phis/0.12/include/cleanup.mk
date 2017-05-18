clean:
	-rm -rf $(cleanfiles)

realclean: clean
	-rm -rf $(cleanfiles) $(realcleanfiles)

distclean:
	-rm -rf $(cleanfiles) $(realcleanfiles) $(distcleanfiles)

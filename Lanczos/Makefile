#CPP = icpc
#CPP_FLAGS = -O3 -axS -ipo
#CPP_FLAGS = -g


#INC_DIR = -I/cvos/shared/apps/acml/4.3.0/ifort64/include  



all:$(LANCZOS_OBJS)

%.o: %.cpp %.h
	$(CPP) $(CPPFLAGS) $(INC_DIR) -c $<

%.o : %.f90 
	$(FC) $(FFLAGS) -c $<
%.o : %.f
	$(FC) $(FFLAGS) -c $<


.PHONY : clean all
clean:
	-rm *.o *~

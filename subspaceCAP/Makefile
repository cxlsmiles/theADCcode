#CPP = icpc

all:$(SUBSPACECAP_OBJS)

%.o: %.cpp %.hpp 
	$(CPP)  $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.c 
	$(CC)  $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

.PHONY : clean all
clean:
	-rm *.o *.mod *~

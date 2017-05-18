

all:$(SCFDATA_OBJS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.c
	$(CC)  $(CFLAGS) $(INC_DIR) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<



.PHONY : clean all
clean:
	-rm *.o *~

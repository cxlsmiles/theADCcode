#CPP = icpc


all:$(ADC2_POL_OBJS)

%.o: %.cpp %.hpp 
	$(CPP)  $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

.PHONY : clean all
clean:
	-rm *.o *.mod *~

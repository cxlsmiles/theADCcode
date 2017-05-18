
#CPP = icpc


#OBJS = adc2_cap_matrix.o zeroth_order.o adc2_cap_analyzer.o 

all:$(ADC2_PROP_OBJS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<


.PHONY : clean all
clean:
	-rm *.o *~



#OBJS = adc2_dip_blocks.o  adc2_matrix.o singlet.o  triplet.o adc2_dip_analyzer.o

all:$(ADC2_OBJS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(OPENMP) $(INC_DIR) -c $<

.PHONY : clean all
clean:
	-rm *.o *~

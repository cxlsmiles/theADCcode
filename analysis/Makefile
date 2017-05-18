

#OBJS = subspacecap_analyzer.o adc2_dip_analyzer.o    adc_analyzer.o    full_cap_analyzer.o    isr_dipole_analyzer.o    

all:$(ANALYSIS_OBJS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(INC_DIR) -c $<

.PHONY : clean all
clean:
	-rm *.o *~

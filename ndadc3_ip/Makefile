#CPP = icpc


all:$(NDADC3_IP_OBJS)

%.o: %.cpp %.hpp 
	$(CPP)  $(CPPFLAGS) $(OPENMP) $(INC_DIR) -c $<

%.o: %.cpp 
	$(CPP)  $(CPPFLAGS) $(OPENMP) $(INC_DIR) -c $<

%.o: %.c 
	$(CC) $(CFLAGS)  -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY : clean all
clean:
	-rm *.o *~

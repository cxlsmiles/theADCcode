#CPP = icpc



all:$(NDADC3_PROP_OBJS)

%.o: %.cpp %.hpp 
	$(CPP) $(CPPFLAGS) $(INC_DIR) -c $<

%.o: %.c 
	$(CC) $(CFLAGS) $(INC_DIR)  -c $<
.PHONY : clean all
clean:
	-rm *.o *~

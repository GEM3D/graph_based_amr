GC= g++ -g

result: main.o compact.o hdf5.o typedefs.h
	$(GC) main.o compact.o hdf5.o -L/usr/local/hdf5_18/1.8.11/lib/ -I/usr/local/hdf5_18/1.8.11/include/ -lhdf5 -o result

main.o: main.cpp  
	$(GC) -c main.cpp
	
hdf5.o: hdf5.c	
	$(GC)  -L/usr/local/hdf5_18/1.8.11/lib/ -I/usr/local/hdf5_18/1.8.11/include/ -lhdf5 -c hdf5.c 
    
compact.o: compact.cpp
	$(GC) -c compact.cpp  -lm  

clean:
	rm *.o *~ result

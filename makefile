CC=g++
CXXFLAGS= -g -v -Wall  $(shell root-config --libs --cflags)
DEPS = plothelper.h  RecoOpAng1.h  SEAobject.h  SEAview_lite_epem.h  SetRealAspectRatio2.h
OBJ = plothelper.o  RecoOpAng1.o  RunSEAview.o  SEAview_lite_epem.o  SetRealAspectRatio2.o 


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS)

RunSEA: $(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS)

clean:
	-rm *.o $(objects) RunSEA


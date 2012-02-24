CXX = mpic++
CXXFLAGS = -g -O2 -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE

all: dc sa_check

dc: dc.o tuple.o
	$(CXX) $(CXXFLAGS) -o dc $^

sa_check: sa_check.o tuple.o
	$(CXX) $(CXXFLAGS) -o sa_check $^

clean:
	rm -f *.o *~

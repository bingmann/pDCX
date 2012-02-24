CC = icc 
MPICC = mpiCC
WARN = 
#-Wall 
#-Winline 
CFLAGS = -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE
# $(WARN)
INCLUDE_DIR =-I /usr/mpich2/include/
LIB_DIR = -L /usr/mpich2/lib/
LIBS =
#-ldmpi
#-lmpich

all: dc sa_check 
	
dc: dc.cpp tuple.o
	$(MPICC) $(CFLAGS) -o dc dc.cpp tuple.o $(WARN)$(LIBS)
sa_check: sa_check.cpp tuple.o
	$(MPICC) $(CFLAGS) -o sa_check sa_check.cpp tuple.o
alltoall: alltoall.cpp tuple.o
	$(MPICC) $(CFLAGS) -o alltoall alltoall.cpp tuple.o
alltoallblank: alltoallblank.cpp tuple.o
	$(MPICC) $(CFLAGS) -o alltoallblank alltoallblank.cpp tuple.o
tuple.o: tuple.cpp
	$(CC) $(CFLAGS) -o tuple.o -c tuple.cpp
merse: merse.c
	$(CC) $(CFLAGS) -o merse merse.c

clean:
	rm -f *.o *~

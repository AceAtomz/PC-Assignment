E1 = scan
E2 = scan_omp
E3 = scan_mpi

C1 = gcc
C2 = mpicc
OMPFLAG = -fopenmp
CFLAGS = -g -Wall
MPIFLAGS = -lm

all: $(E1) $(E2) $(E3)

$(E1): $(E1).c
	$(C1) $(CFLAGS) $(OMPFLAG) $(E1).c -o $(E1)

$(E2): $(E2).c
	$(C1) $(CFLAGS) $(OMPFLAG) $(E2).c -o $(E2)

$(E3): $(E3).c
	$(C2) $(CFLAGS) -o $(E3) $(E3).c $(MPIFLAGS)

clean:
	rm -vf $(E1) $(E2) $(E3)


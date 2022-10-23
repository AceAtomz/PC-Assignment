P1 = scan
P2 = scan_omp
P3 = scan_mpi

E1 = scan
E2 = scan_omp
E3 = scan_mpi

C1 = gcc
C2 = mpicc
OMPFLAG = -fopenmp
CFLAGS = -g -Wall

all: $(E1) $(E2) $(E3)

$(E1): $(P1).c
	$(C1) $(CFLAGS) $(OMPFLAG) $(P1).c -o $(E1)

$(E2): $(P2).c
	$(C1) $(CFLAGS) $(OMPFLAG) $(P2).c -o $(E2)

$(E3): $(P3).c
	$(C2) $(P3).c -o $(E3)

clean:
	rm -vf $(E1) $(E2) $(E3)


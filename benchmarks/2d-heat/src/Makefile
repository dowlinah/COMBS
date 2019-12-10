CFLAGS += -lm -pg
MPICC = mpicc
OMP = -fopenmp

all: serial openmp mpi

serial: serial-2dheat.o
	$(CC) serial-2dheat.o -o serial-2dheat.x $(CFLAGS)

serial-2dheat.o: serial-2dheat.c 
	$(CC) -c serial-2dheat.c $(CFLAGS)

openmp: openmp-2dheat.o
	$(CC) $(OMP) openmp-2dheat.o -o openmp-2dheat.x $(CFLAGS)

openmp-2dheat.o: openmp-2dheat.c 
	$(CC) -c openmp-2dheat.c $(CFLAGS)

mpi: mpi-2dheat.o
	$(MPICC) mpi-2dheat.o -o mpi-2dheat.x $(CFLAGS)

mpi-2dheat.o: mpi-2dheat.c 
	$(MPICC) -c mpi-2dheat.c $(CFLAGS)

bench: serial-2dheat.o
	./bench.sh

clean:
	rm -rf *.o *.x callgrind* time.txt mem.txt massif.out ../callgrind.out* ../time.txt ../mem.txt

PROG=matrix-mpi

all: $(PROG)


$(PROG): main.c
	mpicc -o $(PROG) main.c

bench: $(PROG) 
	./bench.sh ./$(PROG)

clean:
	rm -f $(PROG) callgrind.out* time.txt

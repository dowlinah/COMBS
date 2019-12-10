PROG=ksbenchmark

all: $(PROG)


$(PROG): $(PROG).c
	clang -O3 -o $(PROG) $(PROG).c -L/usr/lib -lfftw3 -lm

bench: $(PROG) 
	./bench.sh ./$(PROG) 50000

clean:
	rm -f $(PROG) *.plt callgrind.out* time.txt massif.out mem.txt
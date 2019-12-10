PROG=RRR

all: $(PROG)


$(PROG): RRR.c
	gcc -O2 RRR.c -lm -L/usr/lib -lfftw3 -o $(PROG) 

bench: $(PROG) 
	./bench.sh

clean:
	rm -f $(PROG) callgrind.out* time.txt results100E sol massif.out mem.txt

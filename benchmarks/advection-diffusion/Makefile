PROG=heat

all: $(PROG)


$(PROG): $(PROG).c
	clang -O3 -o $(PROG) $(PROG).c -lm

bench: $(PROG) 
	./bench.sh ./$(PROG)

clean:
	rm -f $(PROG) *.plt callgrind.out* time.txt massif.out mem.txt
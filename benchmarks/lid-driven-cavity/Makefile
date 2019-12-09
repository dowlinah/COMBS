PROG=Lid-driven

all: $(PROG)


$(PROG): Lid-driven.c
	gcc -o $(PROG) Lid-driven.c

bench: $(PROG) 
	./bench.sh ./$(PROG)

clean:
	rm -f $(PROG) *.plt callgrind.out* time.txt massif.out mem.txt

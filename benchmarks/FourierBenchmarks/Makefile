all:
	$(MAKE) -C FFTW

bench: all
	./bench.sh

clean:
	$(MAKE) clean -C FFTW
	rm -f time.txt callgrind.out* massif.out mem.txt

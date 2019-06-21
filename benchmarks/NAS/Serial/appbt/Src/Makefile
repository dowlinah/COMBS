# Makefile for appbt application

CC = cc
CFLAGS =
LIBS = -lm

#Given i/o sizes are 12 and 64
INPUT_SIZE = 12

APPBT = appbt.o jac.o btrid.o rhs.o

all : 	appbt

appbt :	$(APPBT) appbt.input appbt.output
	$(CC) -o appbt $(CFLAGS) $(APPBT) $(LIBS)

appbt.o :	appbt.c appbt.h

jac.o :		jac.c appbt.h

btrid.o :	btrid.c appbt.h

rhs.o :		rhs.c appbt.h

appbt.input:	appbt.inp$(INPUT_SIZE)
		cp appbt.inp$(INPUT_SIZE) appbt.input

appbt.output:	appbt.out$(INPUT_SIZE)
		cp appbt.out$(INPUT_SIZE) appbt.output

clean:		
		rm -f *.o appbt.input appbt.output appbt
		
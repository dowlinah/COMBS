#
# NAS CFD code: APPSP
#
# Fortran version:
#
# Author: Sisira Weeratunga
#         NASA Ames Research Center
#         (10/25/90)
#
# This is the C version of the Fortran code developed at NASA Ames Research
# Center.  We converted the code to C as a part of our class project for 
# CS 838-3/ChE 562  offered in Spring 1993 by Mark D. Hill, Sangtae Kim and
# Mary Vernon at the University of Wisconsin at Madison.  CS 838-3/ChE 562 
# was an experimental course that brought computer scientists and computation
# scientists together to promote interdisciplinary research.
#
#					June 18, 1993.
#
#					Shubhendu S. Mukherjee (Shubu)
#					shubu@cs.wisc.edu
#					Computer Sciences Department
#
#					Iasonas Moustakis
#					iasonas@luther.che.wisc.edu
#					Chemical Engineering Department
#					
#					University of Wisconsin at Madison
#

SIZE 	= tiny

header 	= common.h input.h

CC	= gcc	
CFLAGS 	= -O2 

TARGET = appsp

OBJS    = adi.o erhs.o error.o exact.o invr.o jac.o l2norm.o main.o pintgr.o \
	  setbv.o setiv.o spenta.o rhs.o tinv.o verify.o 

.SUFFIXES:
.SUFFIXES: .o .c
.c.o: ; $(CC) $(CFLAGS) -c $*.c 

$(TARGET): input.h $(OBJS) 
	$(CC)  $(CFLAGS) -o $(TARGET).$(SIZE) $(OBJS) -lm

input.h: input.$(SIZE).h
	/bin/rm -f input.h
	cp input.$(SIZE).h input.h

install:  $(TARGET).$(SIZE)
	cp $(TARGET).$(SIZE) ../Bin/.

clean:
	/bin/rm -f *.o $(TARGET) input.h \
	$(TARGET).tiny $(TARGET).medium $(TARGET).large 


adi.o:   	adi.c	 ${header}	
erhs.o:	 	erhs.c	 ${header}
error.o: 	error.c  ${header}
exact.o: 	exact.c  ${header}
invr.o:  	invr.c   ${header}
jac.o:	 	jac.c    ${header}
l2norm.o: 	l2norm.c ${header}
main.o:  	main.c   ${header}
pintgr.o:	pintgr.c ${header}
setbv.o: 	setbv.c  ${header}
setiv.o: 	setiv.c  ${header}
spenta.o: 	spenta.c ${header}
rhs.o:	 	rhs.c    ${header}
tinv.o:	  	tinv.c   ${header}
verify.o: 	verify.c ${header}


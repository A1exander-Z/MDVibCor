# Makefile for MDVibCor

CPP = g++

CFLAGS=-c -O2 -Werror

LDFLAGS=-lboost_system -lboost_thread -lgsl -lgslcblas

PROGRAM = mdvibcor

OBJECTS = $(PROGRAM).o atom.o geometry.o pair.o simulation.o

.SUFFIXES: .o .cpp

.cpp.o :
	$(CPP) $(CFLAGS) -o $@ $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CPP) -o $(PROGRAM) $(OBJECTS) $(LDFLAGS)

clean: 
	rm -f *.o $(PROGRAM)


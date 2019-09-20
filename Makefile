
bichsel: bichsel.f90
	gfortran -O2 -o bichsel bichsel.f90

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

simdelta: simdelta.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o simdelta simdelta.cc $(ROOTLIBS)
	@echo "done: simdelta"

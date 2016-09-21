GOPT=-ggdb

all: calcDihedral

calcDihedral: vector_class.o lib.o dihedral.o changeDih.o
	gfortran $(GOPT) -o $@ $^

vector_class.o: vector_class.f90        makefile
	gfortran -c $(GOPT) -o $@ $<

lib.o: lib.f90 vector_class.o  makefile
	gfortran -c $(GOPT) -o $@ $<

dihedral.o: dihedral.f90 vector_class.o  makefile
	gfortran -c $(GOPT) -o $@ $<

changeDih.o: changeDih.f90  dihedral.o lib.o vector_class.o makefile
	gfortran -c $(GOPT) -o $@ $<
        
clean:
	rm -f *.mod *.o $@

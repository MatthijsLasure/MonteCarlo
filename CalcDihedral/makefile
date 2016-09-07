GOPT=-ggdb

all: vector_class.o lib.o dihedral.o changeDih.o
	gfortran $(GOPT) -o $@ $^

vector_class.o: vector_class.f90        makefile
	gfortran -c $(GOPT) -o $@ $<

lib.o: lib.f90  makefile
	gfortran -c $(GOPT) -o $@ $<

dihedral.o: dihedral.f90        makefile
	gfortran -c $(GOPT) -o $@ $<

changeDih.o: changeDih.f90  makefile
	gfortran -c $(GOPT) -o $@ $<
        
clean:
	rm -f *.mod *.o $@

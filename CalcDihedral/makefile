GOPT=-ggdb -fopenmp

calcDihedral: vector_class.o lib.o dihedral.o calcdih.o
	gfortran $(GOPT) -o $@ $^

vector_class.o: vector_class.f90        makefile
	gfortran -c $(GOPT) -o $@ $<

lib.o: lib.f90  makefile
	gfortran -c $(GOPT) -o $@ $<

dihedral.o: dihedral.f90        makefile
	gfortran -c $(GOPT) -o $@ $<

calcdih.o: calcdih.f90  makefile
	gfortran -c $(GOPT) -o $@ $<
        
clean:
	rm -f *.mod *.o $@

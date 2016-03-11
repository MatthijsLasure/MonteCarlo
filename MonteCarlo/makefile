MonteCarlo : randgen.o vector_class.o readConfig.o lib.o LennardJones.o Gaussian.o MonteCarlo.o
	gfortran -ggdb -o MonteCarlo tmp/*.o

MonteCarlo.o : MonteCarlo.f08
	gfortran -ggdb -c -o tmp/MonteCarlo.o MonteCarlo.f08 

LennardJones.o : LennardJones.f08
	gfortran -ggdb -c -o tmp/LennardJones.o LennardJones.f08 

Gaussian.o : Gaussian.f08
	gfortran -ggdb -c -o tmp/Gaussian.o Gaussian.f08 

lib.o : lib.f08
	gfortran -ggdb -c -o tmp/lib.o lib.f08

readConfig.o : readConfig.f08
	gfortran -ggdb -c -o tmp/readConfig.o readConfig.f08

vector_class.o : vector_class.f08
	gfortran -ggdb -c -o tmp/vector_class.o vector_class.f08

randgen.o : randgen.f
	gfortran -ggdb -c -o tmp/randgen.o randgen.f

clean : 
	rm tmp/*.o *.mod MonteCarlo

OpenMP : randgenMP vector_classMP readConfigMP libMP LennardJonesMP GaussianMP MonteCarloMP
	gfortran -ggdb -fopenmp -o MonteCarloOpenMP tmp/*.o

MonteCarloMP : MonteCarlo.f08
	gfortran -ggdb -fopenmp -c -o tmp/MonteCarlo.o MonteCarlo.f08 

LennardJonesMP : LennardJones.f08
	gfortran -ggdb -fopenmp -c -o tmp/LennardJones.o LennardJones.f08 

GaussianMP : Gaussian.f08
	gfortran -ggdb -fopenmp -c -o tmp/Gaussian.o Gaussian.f08 

libMP : lib.f08
	gfortran -ggdb -fopenmp -c -o tmp/lib.o lib.f08

readConfigMP : readConfig.f08
	gfortran -ggdb -fopenmp -c -o tmp/readConfig.o readConfig.f08

vector_classMP : vector_class.f08
	gfortran -ggdb -fopenmp -c -o tmp/vector_class.o vector_class.f08

randgenMP : randgen.f
	gfortran -ggdb -fopenmp -c -o tmp/randgen.o randgen.f

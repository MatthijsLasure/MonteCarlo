MonteCarlo : randgen.o vector_class.o dihedral.o readConfig.o lib.o solmod.o LennardJones.o Gaussian.o MonteCarlo.o
	gfortran -ggdb -fcheck=bounds -o MonteCarlo tmp/*.o

MonteCarlo.o : MonteCarlo.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/MonteCarlo.o MonteCarlo.f08 

LennardJones.o : LennardJones.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/LennardJones.o LennardJones.f08 

Gaussian.o : Gaussian.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/Gaussian.o Gaussian.f08 

lib.o : lib.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/lib.o lib.f08

readConfig.o : readConfig.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/readConfig.o readConfig.f08

vector_class.o : vector_class.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/vector_class.o vector_class.f08

randgen.o : randgen.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/randgen.o randgen.f08

solmod.o : solmod.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/solmod.o solmod.f08

dihedral.o : dihedral.f08
	gfortran -ggdb -c -fcheck=bounds -o tmp/dihedral.o dihedral.f08

clean : 
	rm tmp/*.o *.mod MonteCarlo

OpenMP : randgenMP vector_classMP dihedralMP readConfigMP libMP solmodMP LennardJonesMP GaussianMP MonteCarloMP solmodMP
	gfortran -ggdb -fopenmp -fcheck=bounds -o MonteCarloOpenMP tmp/*.o

MonteCarloMP : MonteCarlo.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/MonteCarlo.o MonteCarlo.f08 

LennardJonesMP : LennardJones.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/LennardJones.o LennardJones.f08 

GaussianMP : Gaussian.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/Gaussian.o Gaussian.f08 

libMP : lib.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/lib.o lib.f08

readConfigMP : readConfig.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/readConfig.o readConfig.f08

vector_classMP : vector_class.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/vector_class.o vector_class.f08

randgenMP : randgen.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/randgen.o randgen.f08

solmodMP : solmod.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/solmod.o solmod.f08

dihedralMP : dihedral.f08
	gfortran -ggdb -fopenmp -c -fcheck=bounds -o tmp/dihedral.o dihedral.f08

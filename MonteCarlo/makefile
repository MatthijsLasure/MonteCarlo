intel : MonteCarloOpenMP
GNU : MonteCarloGNU

clean : 
	-rm tmp-gnu/* MonteCarloGNU MonteCarloOpenMP

turing :
	-mkdir -p turing
	cp MonteCarloOpenMP turing

hopper :
	-mkdir -p hopper
	cp MonteCarloOpenMP hopper

#
# Intel compile
#
ITMP=tmp-intel
IOPT=-O2 -qopenmp -DUSE_PIPES -module $(ITMP)
#IOPT=-O0 -g -qopenmp -DUSE_PIPES -module $(ITMP)
#IOPT=-O0 -g -qopenmp -DDEBUG -DUSE_PIPES -module $(ITMP)
#IOPT=-O0 -g -qopenmp -DDEBUG -module $(ITMP)

MonteCarloOpenMP : $(ITMP)/MCconstants.o $(ITMP)/randgen.o      $(ITMP)/vector_class.o \
                   $(ITMP)/dihedral.o    $(ITMP)/readConfig.o   $(ITMP)/lib.o \
                   $(ITMP)/solmod.o      $(ITMP)/LennardJones.o $(ITMP)/Gaussian.o \
                   $(ITMP)/MonteCarlo.o  $(ITMP)/solmod.o
	ifort $(IOPT) -o $@ $^

$(ITMP)/MCconstants.o  : MCconstants.f90        makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/MonteCarlo.o   : MonteCarlo.F90         $(ITMP)/MCconstants.o  $(ITMP)/vector_class.o \
                         $(ITMP)/LennardJones.o $(ITMP)/Gaussian.o     $(ITMP)/lib.o \
                         $(ITMP)/randgen.o      $(ITMP)/readConfig.o   $(ITMP)/solmod.o       makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/LennardJones.o : LennardJones.f90       makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/Gaussian.o     : Gaussian.F90           $(ITMP)/MCconstants.o  makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/lib.o          : lib.f90                $(ITMP)/MCconstants.o  makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/readConfig.o   : readConfig.f90         $(ITMP)/MCconstants.o  makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/vector_class.o : vector_class.f90       makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/randgen.o      : randgen.f90            makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/solmod.o       : solmod.F90             $(ITMP)/MCconstants.o  $(ITMP)/randgen.o \
                         $(ITMP)/dihedral.o     $(ITMP)/lib.o          $(ITMP)/vector_class.o  makefile
	ifort -c $(IOPT) -o $@ $< 

$(ITMP)/dihedral.o     : dihedral.f90           $(ITMP)/MCconstants.o  $(ITMP)/vector_class.o  makefile
	ifort -c $(IOPT) -o $@ $< 


#
# gfortran compile with OpenMP
#
GTMP=tmp-gnu
GOPT= -fopenmp -DUSE_PIPES -J$(GTMP)
MonteCarloGNU : $(GTMP)/MCconstants.o $(GTMP)/randgen.o      $(GTMP)/vector_class.o \
                $(GTMP)/dihedral.o    $(GTMP)/readConfig.o   $(GTMP)/lib.o \
                $(GTMP)/solmod.o      $(GTMP)/LennardJones.o $(GTMP)/Gaussian.o \
                $(GTMP)/MonteCarlo.o  $(GTMP)/solmod.o
	gfortran $(GOPT) -o $@ $^

$(GTMP)/MCconstants.o  : MCconstants.f90        makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/MonteCarlo.o   : MonteCarlo.F90         $(GTMP)/MCconstants.o  $(GTMP)/vector_class.o \
                         $(GTMP)/LennardJones.o $(GTMP)/Gaussian.o     $(GTMP)/lib.o \
                         $(GTMP)/randgen.o      $(GTMP)/readConfig.o   $(GTMP)/solmod.o       makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/LennardJones.o : LennardJones.f90       makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/Gaussian.o     : Gaussian.F90           $(GTMP)/MCconstants.o  makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/lib.o          : lib.f90                $(GTMP)/MCconstants.o  makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/readConfig.o   : readConfig.f90         $(GTMP)/MCconstants.o  makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/vector_class.o : vector_class.f90       makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/randgen.o      : randgen.f90            makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/solmod.o       : solmod.f90             $(GTMP)/MCconstants.o  $(GTMP)/randgen.o \
                         $(GTMP)/dihedral.o     $(GTMP)/lib.o          $(GTMP)/vector_class.o  makefile
	gfortran -c $(GOPT) -o $@ $< 

$(GTMP)/dihedral.o     : dihedral.f90           $(GTMP)/MCconstants.o  $(GTMP)/vector_class.o  makefile
	gfortran -c $(GOPT) -o $@ $< 



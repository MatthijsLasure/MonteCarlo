all: GNU

LIBS+=MC_lib.so
MODS+=/home/organica/matthijs/libsrc/
OPS+=DUSE_PIPES

F90_SRCS += \
Gaussian.F90 \
LennardJones.f90 \
MonteCarlo.F90 \
readConfig.f90 \
solmod.F90

OBJS += \
Gaussian.o \
LennardJones.o \
MonteCarlo.o \
readConfig.o \
solmod.o

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f90
        @echo "Building file: $<"
        @echo "Invoking: GNU Fortran Compiler"
        gfortran -fopenmp -I$(MODS) -Wall -c -o "$@" "$<" $(LIBS)
        @echo "Finished building: $<"
        @echo " "

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.F90
        @echo "Building file: $<"
        @echo "Invoking: GNU Fortran Compiler"
        gfortran -fopenmp -I$(MODS) -D$(OPS) -Wall -c -o "$@" "$<" $(LIBS)
        @echo "Finished building: $<"
        @echo " "

GNU: $(OBJS)
        @echo "Building target: $@"
        @echo "Invoking: GNU Fortran Linker"
        gfortran -fopenmp -Wall -I$(MODS) -o "MonteCarloGNU" $(OBJS) $(USER_OBJS) $(LIBS)
        @echo "Finished building target: $@"
        @echo " "

# Objects
Gaussian.o: Gaussian.F90
LennardJones.o: LennardJones.f90
MonteCarlo.o: MonteCarlo.F90 Gaussian.o LennardJones.o readConfig.o solmod.o
readConfig.o: readConfig.f90
solmod.o: solmod.F90

clean:
        -$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) MonteCarloGNU
        -@echo " "

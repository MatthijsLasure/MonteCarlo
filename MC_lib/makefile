all: library

F90_SRCS += \
dihedral.f90 \
lib.f90 \
lib_read.f90 \
MCconstants.f90 \
randgen.f90 \
vector_class.f90

OBJS += \
dihedral.o \
lib.o \
lib_read.o \
MCconstants.o \
randgen.o \
vector_class.o


%.o: %.f90
	@echo "Building file: $<"
	@echo "Invoking: GNU Fortran Compiler"
	gfortran -c -fPIC -o "$@" "$<"
	@echo "Finished building: $<"
	@echo " "


library: $(OBJS)
	@echo "Building target: $@"
	@echo "Invoking: GNU Fortran Linker"
	gfortran -shared -o "MC_lib.so" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo "Finished building target: $@"
	@echo " "

vector_class.o: vector_class.f90

randgen.o: randgen.f90

MCconstants.o: MCconstants.f90

lib.o: lib.f90 vector_class.o randgen.o

lib_read.o: lib_read.f90 vector_class.o

dihedral.o: dihedral.f90 MCconstants.o vector_class.o

clean:
	-$(RM) $(OBJS) MC_lib.dll
	-@echo " "

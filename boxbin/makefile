all: boxbin


boxbin: $(OBJS)
	@echo "Building target: $@"
	@echo "Invoking: GNU Fortran Linker"
	gfortran  -o "boxbin" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo "Finished building target: $@"
	@echo " "

	# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) boxbin
	-@echo " "

.PHONY: all clean dependents
.SECONDARY:

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../boxbin.f90 \
../lib.f90 \
../vector_class.f90 

OBJS += \
./boxbin.o \
./lib.o \
./vector_class.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo "Building file: $<"
	@echo "Invoking: GNU Fortran Compiler"
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo "Finished building: $<"
	@echo " "

boxbin.o: ../boxbin.f90 lib.o vector_class.o

lib.o: ../lib.f90 vector_class.o

vector_class.o: ../vector_class.f90

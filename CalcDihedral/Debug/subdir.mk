################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../calcdih.f90 \
../dihedral.f90 \
../vector_class.f90 

F08_SRCS += \
../lib.f08 

OBJS += \
./calcdih.o \
./dihedral.o \
./lib.o \
./vector_class.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

calcdih.o: ../calcdih.f90 dihedral.o lib.o vector_class.o

dihedral.o: ../dihedral.f90 vector_class.o

%.o: ../%.f08
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

lib.o: ../lib.f08 vector_class.o

vector_class.o: ../vector_class.f90



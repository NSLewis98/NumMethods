FC=gfortran
FFLAGS=-g
LINK= -L ../lapack-3.9.0 -llapack -lrefblas -lblas
SRC= ../src/ModuleAngularMomentum.f90
DEBUG_F= -g -fbacktrace -Wall -fcheck=all
MAIN=modIntegrate.f90

build:
	${FC} ${SRC}  ${MAIN} -w -o main 

debug:	
	${FC} ${FFLAGS} ${DEBUG_F} ${SRC} ${LINK}
	
clean:
	@rm -f *.dat *.mod *.o *.out main

buildMod:
	${FC} ${SRC} -c ${MODULE} ${LINK}

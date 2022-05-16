FC=gfortran
FFLAGS=-g
LINK= -L ../lapack-3.9.0 -llapack -lrefblas -lblas
SRC= ../src/ModuleAngularMomentum.f90 ../src/modSphericalIntegral.f90
DEBUG_F= -g -fbacktrace -Wall -fcheck=all
MAIN=ProgLebQuad.f90

build:
	${FC} ${SRC}  ${MAIN} -w -o main 

debug:	
	${FC} ${FFLAGS} ${DEBUG_F} ${SRC} ${LINK}
	
clean:
	@rm -f *.dat *.mod *.o *.out 	

buildMod:
	${FC} ${SRC} -c ${MODULE} ${LINK}

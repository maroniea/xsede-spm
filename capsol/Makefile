FC=gfortran
FFLAGS=-llapack -g


%.o: %.f90
	$(FC) -c $(FFLAGS) $<


test: test_capsolpy.o CapSolPy.o
	$(FC) $(FFLAGS) test_capsolpy.o CapSolPy.o -o test
# load python on OSC 
module load python/3.7-2019.10
#create python/fortran file
f2py --f90flags=-mkl --fcompiler=intel -c capsolpy.pyf CapSolPy.f90

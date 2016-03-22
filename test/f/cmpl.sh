gfortran -fopenmp -g -o test test.f90 ../../fort/sat_f.a ../../sat/sat.a ../../gp/gp.a ../../Lbfgsb.3.0/liblbfgsb.a ~/lapack-3.5.0/liblapack.a ~/lapack-3.5.0/BLAS/SRC/librefblas.a -lm

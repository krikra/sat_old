gcc -g -O0 -o gp ~/lapack-3.5.0/BLAS/SRC/librefblas.a ~/lapack-3.5.0/liblapack.a sample.c gp.a -I ../include/ -lm -llapack -lblas

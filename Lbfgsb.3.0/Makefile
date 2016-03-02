FC = gfortran

#FFLAGS = -O3 -Wall -fbounds-check -g -Wno-uninitialized 
FFLAGS = -O3 -g

DRIVER1_77 = driver1.f
DRIVER2_77 = driver2.f
DRIVER3_77 = driver3.f

DRIVER1_90 = driver1.f90
DRIVER2_90 = driver2.f90
DRIVER3_90 = driver3.f90

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f

OBJ = lbfgsb.o linpack.o timer.o

lib: liblbfgsb.a

all :  lbfgsb_77_1 lbfgsb_77_2 lbfgsb_77_3 lbfgsb_90_1 lbfgsb_90_2 lbfgsb_90_3 


liblbfgsb.a: $(OBJ)
	ar rcs $@ $^

%.o : %.f
	$(FC) $(FFLAGS) -c $^

lbfgsb_77_1 : $(DRIVER1_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_1

lbfgsb_77_2 : $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_2

lbfgsb_77_3 : $(DRIVER3_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_3

lbfgsb_90_1 : $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_1

lbfgsb_90_2 : $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_2

lbfgsb_90_3 : $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_3

clean:
	rm -f *.o
	rm -f *.a

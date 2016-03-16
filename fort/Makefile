include ../Makefile.in

#INC=../include
#CFLAGS= -O0 -lm -fopenmp -g

OBJ=sat_f_interface.o sat_f_interface.o


sat_f.a: ${OBJ}
	ar rcs $@ $^

${OBJ}:

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

%.o: %.f90
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o
	rm -f *.mod

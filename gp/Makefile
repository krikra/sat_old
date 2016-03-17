include ../Makefile.in


OBJ=gp_kernel.o gp_power.o gp_misc.o gp_mem.o

gp.a: ${OBJ}
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o

include ../Makefile.in

OBJ=sat_sys.o sat_mem.o ippe_id.o

sat.a: ${OBJ}
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o

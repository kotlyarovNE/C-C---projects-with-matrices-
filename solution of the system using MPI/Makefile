CFLAGS	= -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wcast-qual -Wno-suggest-attribute=format


all :  a.out

a.out : main.o input.o otlad.o output.o nevazka.o functions.h 
	mpicc main.o input.o otlad.o  output.o nevazka.o functions.h -lm -o ./a.out   -O3

main.o : main.c functions.h
	mpicc -c $(CFLAGS) main.c   		-O3
	
input.o : input.c functions.h
	mpicc -c $(CFLAGS) input.c 		-O3

otlad.o : otlad.c functions.h
	mpicc -c $(CFLAGS) otlad.c 		-O3

output.o : output.c functions.h
	mpicc -c $(CFLAGS) output.c 		-O3

nevazka.o : nevazka.c functions.h
	mpicc -c $(CFLAGS) nevazka.c 		-O3

clean:
	rm -rf *.o ./a.out

CC = g++
CFLAGS = -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

all: a.out

a.out: main.o f.o input.o output.o solution.o function.h
	$(CC) $(CFLAGS) main.o f.o input.o output.o solution.o function.h -o a.out

main.o: main.cpp
	$(CC) -c $(CFLAGS) -c main.cpp

f.o: f.cpp 
	$(CC) -c $(CFLAGS) -c f.cpp

input.o: input.cpp 
	$(CC) -c $(CFLAGS) -c input.cpp

output.o: output.cpp
	$(CC) -c $(CFLAGS) -c output.cpp

solution.o: solution.cpp
	$(CC) -c $(CFLAGS) -c solution.cpp

clean:
	rm -f *.o a.out
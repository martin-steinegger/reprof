# the compiler to use.
CC=gcc
# options I'll pass to the compiler.
CFLAGS=-x c -m64 -std=gnu99 -Wno-trigraphs -Wno-deprecated -Wno-char-subscripts -O3 -Wall -funroll-loops -msse2 -msse3 -lm -fvisibility=hidden -I. -Iinclude

OBJECTS=main.o

all: reprof

reprof: $(OBJECTS)
	$(CC) $(OBJECTS) -static ../fannorg/src/libfann.a -lz -lm -o reprof 

main.o: main.c
	$(CC) $(CFLAGS) -c main.c 

clean:
	rm -rf *.o reprof

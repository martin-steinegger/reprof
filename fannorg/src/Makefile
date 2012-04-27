# the compiler to use.
CC=gcc
# options I'll pass to the compiler.
CFLAGS=-x c -m64 -std=gnu99 -Wno-trigraphs -Wno-char-subscripts -O4 -Wall -funroll-loops -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -lm  -fvisibility=hidden -I. -Iinclude

OBJECTS=floatfann.o doublefann.o fann.o fann_error.o fann_io.o fann_train.o fann_train_data.o fixedfann.o fann_cascade.o

all: fann

fann: $(OBJECTS)
	ar cru libfann.a $(OBJECTS)
	ranlib libfann.a
floatfann.o: floatfann.c
	$(CC) $(CFLAGS) -c floatfann.c 

doublefann.o: doublefann.c
	$(CC) $(CFLAGS) -c doublefann.c
		
fann.o: fann.c
	$(CC) $(CFLAGS) -c fann.c

fann_error.o: fann_error.c
	$(CC) $(CFLAGS) -c fann_error.c
		
fann_io.o: fann_io.c
	$(CC) $(CFLAGS) -c fann_io.c

fann_train.o: fann_train.c
	$(CC) $(CFLAGS) -c fann_train.c
		
fann_train_data.o: fann_train_data.c
	$(CC) $(CFLAGS) -c fann_train_data.c

fixedfann.o: fixedfann.c
	$(CC) $(CFLAGS) -c fixedfann.c

fann_cascade.o: fann_cascade.c
	$(CC) $(CFLAGS) -c fann_cascade.c
clean:
	rm -rf *.o libfann.a
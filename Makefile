CC=gcc
CFLAGS= -O3 -march=native -std=c11
LIBS=-lgsl -lm -lgslcblas

all: noisy

noisy: noisy.c
	$(CC) $(CFLAGS) -o noisy noisy.c $(LIBS)

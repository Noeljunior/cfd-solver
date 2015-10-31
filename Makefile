CFLAGS =  -m64 -Wall
LIBS   = -lm -lgsl -lgslcblas
OPTS   = -march=native -mtune=native -O3 -funroll-loops
SRC=src/*
all:
	gcc $(CFLAGS) $(LIBS) $(OPTS) -I$(SRC) $(SRC) -o cfd

run: all
	./cfd

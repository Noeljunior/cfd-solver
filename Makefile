CFLAGS = -march=native -mtune=native
LIBS   =-lm
SRC=src/*
all:
	gcc $(CFLAGS) $(LIBS) -I$(SRC) $(SRC) -o cfd

run: all
	./cfd

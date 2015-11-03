CFLAGS   = -Wall
LIBS     = -lm -lgsl -lgslcblas -lGL -lSDL2 -lGLEW -lGL -lGLU -lpthread
OPTS     = -march=native -mtune=native -m64 -O3 -funroll-loops

SRC      = src/
OBJ      = obj/
OUT      = cfd

# list of object that will be compiling with optimazion flags
OPT_OBJS = cfdsolver.o meshviewer.o
# and the auto-list of the remaining objects
OBJS = $(filter-out $(OPT_OBJS),$(patsubst %.c,%.o,$(subst $(SRC),,$(wildcard $(SRC)*.c))))


$(OUT): $(addprefix ${OBJ},$(OBJS)) $(addprefix ${OBJ},$(OPT_OBJS))
	gcc $(LIBS) ${OBJ}*.o -o $@


$(addprefix ${OBJ},$(OPT_OBJS)): OPT_OBJS_FLAGS := $(OPTS)
${OBJ}%.o: $(SRC)%.c
	@mkdir -p ${OBJ}
	gcc $(CFLAGS) $(OPT_OBJS_FLAGS) -o $@ -c $<

clean:
	rm -f cfd
	rm -fr ${OBJ}

run: $(OUT)
#	./$(OUT) mesh/naca0012_2261VERT.mesh -p1
	./$(OUT) mesh/naca0012_2261VERT.mesh -O2 -C4.5 -I5 -T0.001

force: clean ${OUT}


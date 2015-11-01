/*
 *  This is an implementation of a CFD solver ussing FVM aproach with a
 *  classifier.
 *
 *  This implementation is VERY BASED on work of
 *  Juan Eduardo Casavilca Silva <casavilca.je@gmail.com>.
 *  Everything inside \/ marks, are direct reference to Juan's work.
 *
 *  Although this implementation might look very different, the credits of
 *  original work go to Juan Casavilca.
 *
 *  (C) Tiago Levita <tlevita@student.dei.uc.pt>
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                          DATA STRUCTURE
 *                                           and typedefs
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  Main data structure, as anonymous
 */
typedef struct cfds_mesh cfds_mesh;

/*
 *  Problem-related and solver-behaviour settings
 */
typedef struct cfds_args {
    int     order;              /* Reconstruction order */
    double  cfl;                /* CFL Condiction */
    int     max_iterations;     /* Maximum iterations */
    double  nr_threashold;      /* residue norm threshold to stop */

    char    verbose;            /* if to show something more */
    char    quiet;              /* show nothing */
    char    showinner;          /* if to show detailed results in each rungekutta iterations */
    char    fclassify;          /* if quiet is set and this is set to, then it outputs only the classification */
} cfds_args;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                       FUNCTION DECLARATION
 *                                       of external interface
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Init and get a pointer to the internal data struct
 */
cfds_mesh *     cfds_init(cfds_args ina, double ** vertices, int sizev, int ** edges, int sizee, int ** triangles, int sizet);

/*
 * Solve the problem CFD problem
 */
void            cfds_solve(cfds_mesh * inm);

/*
 * Free up everything else on memory
 */
void            cfds_free(cfds_mesh * inm);



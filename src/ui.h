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

#define PROGRAM_NAME    "msc-cfd"
#define PROGRAM_VERSION "1.0"
#define BUG_REPORT      "tlevita@student.dei.uc.pt"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                          DATA STRUCTURE
 *                                           and typedefs
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  The arguments parsed
 */
typedef struct ui_args {
    char    *mesh_file;         /* what is this */
    /*char    *parameters_file;
    char    *log_file;*/
    char    verbose;            /* if to show something more */
    char    quiet;              /* show nothing */
    char    showinner;          /* if to show detailed results in each rungekutta iterations */
    char    fclassify;          /* if quiet is set and this is set to, then it outputs only the classification */
    char    rstdin;             /* read from stdin and not from a file */
    char    mftype;             /* mesh-file type; 0: is auto; 1: force new; 2: force legacy */

    int     order;              /* reconstruction order */
    double  cfl;                /* CFL condiction */
    int     max_iterations;     /* max iterations */
    double  nr_threashold;      /* residue norm threshold to stop */

    /*char    f_ds,
            f_plot,
            f_pc,
            f_residue;*/
    int     testp;              /* the number of pre-defined problem related settings */

    int     mandatory;          /* inner var to control the mandatory arguments */

    
} ui_args;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                       FUNCTION DECLARATION
 *                                       of external interface
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Init and get a pointer to the internal data struct
 */
ui_args *       ui_readargs(int argc, char **argv);

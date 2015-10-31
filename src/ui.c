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

#include "ui.h"

#include "infocolour.h"

#include <stdlib.h>
#include <argp.h>
#include <error.h>
#include <limits.h>

/*
 *  This file name, used in INFO/WARN/ERR msgs
 */
/*static char *MOD = "UI";*/

/*                      TODO list

    * check for the list of available options
        are they needed?
    * check for restrictions
    * add ability to read from stdin
    * add ability 
    * separete the colours to the same file of INFO/WARN/ERR

    * * NOT FOR NOW / NOT IMPORTANT/RELEVANT **


    * * ALREADY DONE / PRETTY MUCH DONE **


*/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             FUNCTION DECLARATION and globals definitions
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
static error_t      parse_opt(int key, char *arg, struct argp_state *state);
static int          validate_args(ui_args * args);
void static         print_args(ui_args * args, struct argp argp);

static int          cstrtoi_ordie(char * arg, char key, struct argp_state *state);
static double       cstrtod_ordie(char * arg, char key, struct argp_state *state);

/*
 * argp globals: program version and bug report address
 */
const char  *argp_program_version     = PROGRAM_NAME" v"PROGRAM_VERSION;
const char  *argp_program_bug_address = "<"BUG_REPORT">";

/*
 * description text
 */
//#define UNDR "MERDA"
static char doc[] =
  BOLD PROGRAM_NAME" v"PROGRAM_VERSION RESET" -- an efficient CFD solver\
\vThis solver comes, by now, "UND"without any waranty of anything"RESET".\
 It might work some day, or "BOLD"not"RESET".";
/*
 * usage messages
 */
static char args_doc[] =
    "<mesh_file> -O <order> -C <cfl> -I <max_iterations> -T <threshold>\n"
    "-s -O <order> -C <cfl> -I <max_iterations> -T <threshold>\n"
    "<mesh_file> -p <test-problem-id>\n"
    "-s -p <test-problem-id>";

/*
 * list options
 */
struct argp_option options[] = {
    {0,0,0,0,             "Main options"},
    /*{"parameters",        'p', "FILE", OPTION_NO_USAGE,     "A FILE containing the input parameters of the solver", 0},*/
    /*{"log",               'l', "FILE", OPTION_NO_USAGE,     "Output to FILE instead of standard output"},*/
    {"force-classify",    'f', 0     , OPTION_NO_USAGE,     "If quieted, output the classification"},
    {"force-new",         'n', 0     , OPTION_NO_USAGE,     "Tell parser that the mesh file is the new format"},
    {"force-legacy",      'o', 0     , OPTION_NO_USAGE,     "Tell parser that the mesh file is the legacy format; by default, it auto-detects"},
    {"details",           'd', 0     , OPTION_NO_USAGE,     "Show more details during Runge-Kutta iterations"},
    {"stdin",             's', 0     , OPTION_NO_USAGE,     "Read from standard input instead of a file"},
    {"verbose",           'v', 0     , OPTION_NO_USAGE,     "Produce verbose output"},
    {"quiet",             'q', 0     , OPTION_NO_USAGE,     "Show nothing"},

    {0,0,0,OPTION_DOC,    "Problem related options"},
    {"p-order",           'O', "O"   , '0'            ,     "Reconstruction order"},
    {"p-cfl",             'C', "CFL" , '0'            ,     "CFL Condiction"},
    {"p-iterations",      'I', "I"   , '0'            ,     "Maximum iterations"},
    {"p-threshold",       'T', "T"   , '0'            ,     "Norm residue threshold"},
    {0,0,0,0, ""}, /* non-visual group */
    /*{"p-print-ds",        'D', 0     , '0'            ,     "Print data structure"},
    {"p-print-plot",      'P', 0     , '0'            ,     "Print data to plot"},
    {"p-print-pc",        'A', 0     , '0'            ,     "Print pressure coeficients (Cp)"},
    {"p-print-residue",   'R', 0     , '0'            ,     "Print norm residue"},*/
    {0,0,0,0, ""}, /* non-visual group */
    {"p-test-problem",    'p', "ID"  , '0'            ,     "Use pre-defined values"},
    /* the last one is always like this */
    { 0 }
};

/*
 * max mandatory = 2^(number of mandatory options)
 */
int max_mandatory = 16;



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             EXTERNAL INTERFACE
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
ui_args * ui_readargs(int argc, char **argv) {
    ui_args * args = calloc(1, sizeof(ui_args));
    struct argp argp = { options, parse_opt, args_doc, doc };

    /* loop the arguments and options */
    argp_parse(&argp, argc, argv, 0, 0, args);

    /* validate parsed options */
    validate_args(args);

    /* print them if allowed by verbosity/quiety */
    print_args(args, argp);

    /* return them back to the caller */
    return    args;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                                  UI
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    ui_args *args = state->input;

    switch (key) {
        /*case 'p': args->parameters_file  = arg;                            break;
        case 'l': args->log_file         = arg;                            break;*/
        case 's': args->rstdin           = 1;                              break;
        case 'f': args->fclassify        = 1;                              break;
        case 'd': args->showinner        = 1;                              break;
        case 'n': args->mftype           = 1;                              break;
        case 'o': args->mftype           = 2;                              break;
        case 'q': args->quiet            = 1;                              break;
        case 'v': args->verbose          = 1;                              break;

        case 'O': args->order            = cstrtoi_ordie(arg, key, state);
                  args->mandatory       |= 1;                              break;
        case 'C': args->cfl              = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 2;                              break;
        case 'I': args->max_iterations   = cstrtoi_ordie(arg, key, state);
                  args->mandatory       |= 4;                              break;
        case 'T': args->nr_threashold  = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 8;                              break;
        /*case 'D': args->f_ds             = 1;                              break;
        case 'P': args->f_plot           = 1;                              break;
        case 'A': args->f_pc             = 1;                              break;
        case 'R': args->f_residue        = 1;                              break;*/
        case 'p': args->testp            = cstrtoi_ordie(arg, key, state);
                  args->mandatory       |= max_mandatory;                  break;

        case ARGP_KEY_ARG:
            if (state->arg_num >= (1 - args->rstdin)) { /* Too many arguments. */
                printf("\nToo many arguments\n\n");
                argp_usage(state);
            }
            args->mesh_file = arg;
            break;

        case ARGP_KEY_NO_ARGS:
            /* If we reach this, check if user asked to use stdin as input file */
            if (args->rstdin) {
                /* Ok... */
            } else { /* Not ok! */
                printf("\nToo few arguments -- missing mesh file\n\n");
                argp_usage (state);
            }
            break;

        case ARGP_KEY_END:
            /* Check for mandatory options */
            if (args->mandatory < (max_mandatory - 1)) {
                printf("\nToo few options\n\n");
                argp_usage (state);
            }
            break;

        default:
            return ARGP_ERR_UNKNOWN; /* TODO: FATAL ERROR */
    }
    return 0;
}

static int validate_args(ui_args * args) {
    /* check if mesh file exists */
    FILE *f;
    if (!args->rstdin) {
        if ((f = fopen(args->mesh_file, "r")))
            fclose(f);
        else {
            fprintf(stderr, ERRT"ERROR"RESET": file '%s' could not be opened!\n", args->mesh_file);
            exit(1);
        }
    }
    /* check if log file exists */
    /*if (args->log_file && (f = fopen(args->log_file, "r")))
        fclose(f);
    else if (args->log_file) {
        fprintf(stderr, ERRT"ERROR"RESET": file '%s' could not be opened!\n", args->log_file);
        exit(1);
    }*/

    /* check if parameters file exists */
    /*if (args->parameters_file && (f = fopen(args->parameters_file, "r")))
        fclose(f);
    else if (args->log_file) {
        fprintf(stderr, ERRT"ERROR"RESET" file '%s' could not be opened!\n", args->parameters_file);
        exit(1);
    }*/

    /* Check for a preset */
    switch (args->testp) {
        case 0: /* no test problem */ break;
        default:
            printf(WARNT"WARN"RESET": test suit no recognized -- %d. Using the first one\n", args->testp);
        case 1: /* so as Juan Casavilca has on his first example */
            args->order           = 2;
            args->cfl             = 4.5;
            args->max_iterations  = 25;
            args->nr_threashold   = 1e-10;
            break;
            
    }

    /* TODO
     * check for bounderies
     */


    return 0;
}

void static print_args(ui_args * args, struct argp argp) {
    /* TODO: Print input settings ONLY IF VERBOSE AND QUIET ALLOWS THAT */

    printf("\n");
    argp_help (&argp, stdout, ARGP_HELP_PRE_DOC, 0);
    printf("\n");


    printf(BOLD"Input mesh file:"RESET" '%s'\n", args->mesh_file);
    switch (args->mftype) {
        case 0:     printf("\tIt will detect the file format\n\n");
            break;
        case 1:     printf("\tForcing read the file in the "UND"new format"RESET"\n\n");
            break;
        case 2:     printf("\tForcing read the file in the "UND"legacy format"RESET"\n\n");
            break;
    }


    printf(BOLD"Starting the solver with:"RESET"\n\
    Order             : %d\n\
    CFL               : %e\n\
    Max interations   : %d\n\
    Residual threshold: %e\n\n",
    args->order, args->cfl, args->max_iterations, args->nr_threashold);

    printf(BOLD"Verbose info: "RESET"%s%s%s\b\b \n",
    (/*!args->f_ds && !args->f_plot && !args->f_pc && !args->f_residue && */!args->verbose) ? "normal, " : "",
    args->verbose   ? "more verbose, " : "",
    args->showinner ? "show details during iterations, " : ""/*,
    args->f_ds      ? "data structure, " : "",
    args->f_plot    ? "plot data, " : "",
    args->f_pc      ? "pressure coeficients (Cp), " : "",
    args->f_residue ? "norm residue, " : ""*/);

    printf("\n");
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                              HELPERS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * String to integer or exit
 */
int cstrtoi_ordie(char * arg, char key, struct argp_state *state) {
    errno = 0;
    char * c;
    long val = strtol(arg, &c, 0);
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
            || (errno != 0 && val == 0)
            || ((long) INT_MIN > val || val > (long) INT_MAX)
            || (arg == c)) {
        printf("\n%c -- invalid argument: '%s'\n\n", key, arg);
        argp_usage(state);
    }

    return (int) val;
}

/*
 * String to double or exit
 */
double cstrtod_ordie(char * arg, char key, struct argp_state *state) {
    errno = 0;
    char * c;
    double val = strtod(arg, &c);
    if ((errno != 0) || (arg == c)) {
        printf("\n%c -- invalid argument: '%s'\n\n", key, arg);
        argp_usage(state);
    }

    return val;
}



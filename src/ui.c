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


    * * NOT FOR NOW / NOT IMPORTANT/RELEVANT **
    *** new opts
        number of digits to residue AND/OR objective function


    * * ALREADY DONE / PRETTY MUCH DONE **
    * add ability to read from stdin
    * separete the colours to the same file of INFO/WARN/ERR
    * check for restrictions
    * check for the list of available options


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
    {"force-classify",    'f', 0     , OPTION_NO_USAGE,     "If quieted, outputs the last solution's details"},
    {"force-new",         'n', 0     , OPTION_NO_USAGE,     "Tell parser that the mesh file is the new format"},
    {"force-legacy",      'o', 0     , OPTION_NO_USAGE,     "Tell parser that the mesh file is the legacy format; by default, it auto-detects"},
    {"details",           'd', "N"   , OPTION_NO_USAGE | OPTION_ARG_OPTIONAL, "Show more details during Runge-Kutta iterations; if N is specified, show it each N iterations; N = -1 shows only the solution's detail"},
    {"grapgics",          'g', "N"   , OPTION_NO_USAGE | OPTION_ARG_OPTIONAL, "Use graphic output; if N is specified, show it each N iterations; N = -1 shows only the solution's graphics"},
    {"stdin",             's', 0     , OPTION_NO_USAGE,     "Read from standard input instead of a file"},
    {"verbose",           'v', 0     , OPTION_NO_USAGE,     "Produce verbose output"},
    {"quiet",             'q', 0     , OPTION_NO_USAGE,     "Show nothing"},

    {0,0,0,OPTION_DOC,    "Problem related options"},
    {"p-angle",           'A', "A"   , '0'            ,     "Angle of atack in degrees"},
    {"p-mach",            'M', "M"   , '0'            ,     "Mach speed"},
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
int max_mandatory = 1 << 6;



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
        case 'd': args->showdetails      = arg == NULL ? -1 : cstrtoi_ordie(arg, key, state);   break;
        case 'g': args->showgraphics     = arg == NULL ? -1 : cstrtoi_ordie(arg, key, state);   break;
        case 'n': args->mftype           = 1;                              break;
        case 'o': args->mftype           = 2;                              break;
        case 'q': args->quiet            = 1;                              break;
        case 'v': args->verbose          = 1;                              break;

        case 'A': args->angle            = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 1 << 0;                         break;
        case 'M': args->mach             = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 1 << 1;                         break;
        case 'O': args->order            = cstrtoi_ordie(arg, key, state);
                  args->mandatory       |= 1 << 2;                         break;
        case 'C': args->cfl              = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 1 << 3;                         break;
        case 'I': args->max_iterations   = cstrtoi_ordie(arg, key, state);
                  args->mandatory       |= 1 << 4;                         break;
        case 'T': args->nr_threashold    = cstrtod_ordie(arg, key, state);
                  args->mandatory       |= 1 << 5;                         break;
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
                printf("\nToo few options -- misses: ");
                if (!(args->mandatory & (1 << 0))) printf("-A --p-angle, ");
                if (!(args->mandatory & (1 << 1))) printf("-M --p-mach, ");
                if (!(args->mandatory & (1 << 2))) printf("-O --p-order, ");
                if (!(args->mandatory & (1 << 3))) printf("-C --p-cfl, ");
                if (!(args->mandatory & (1 << 4))) printf("-I --p-iterations, ");
                if (!(args->mandatory & (1 << 5))) printf("-T --p-threshold, ");
                printf("\b\b  \n\n");
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
            printf(WARNT"WARN"RESET": test suit not recognized -- %d. Using the first one\n", args->testp);
        case 1: /* so as Juan Casavilca has on his first example */
            if (!(args->mandatory & (1 << 0))) args->angle           = 1.25;
            if (!(args->mandatory & (1 << 1))) args->mach            = 0.8;
            if (!(args->mandatory & (1 << 2))) args->order           = 2;
            if (!(args->mandatory & (1 << 3))) args->cfl             = 3.5;
            if (!(args->mandatory & (1 << 4))) args->max_iterations  = 1000000;
            if (!(args->mandatory & (1 << 5))) args->nr_threashold   = 1e-10;
            break;
        case 2: /* Another example */
            if (!(args->mandatory & (1 << 0))) args->angle           = 3.00;
            if (!(args->mandatory & (1 << 1))) args->mach            = 0.8;
            if (!(args->mandatory & (1 << 2))) args->order           = 2;
            if (!(args->mandatory & (1 << 3))) args->cfl             = 2.5;
            if (!(args->mandatory & (1 << 4))) args->max_iterations  = 100000;
            if (!(args->mandatory & (1 << 5))) args->nr_threashold   = 1e-6;
            break;


        case 666: /* try it yourself ;) */
            if (!(args->mandatory & (1 << 0))) args->angle           = 45.00;
            if (!(args->mandatory & (1 << 1))) args->mach            = 0.8;
            if (!(args->mandatory & (1 << 2))) args->order           = 2;
            if (!(args->mandatory & (1 << 3))) args->cfl             = 5.0;
            if (!(args->mandatory & (1 << 4))) args->max_iterations  = 100000;
            if (!(args->mandatory & (1 << 5))) args->nr_threashold   = 1e-6;
            break;
    }

    /* check for auto show*'s */
    if (args->showdetails < 0)
        args->showdetails = (args->max_iterations / 25) < 1 ? 1 : (args->max_iterations / 25);
    if (args->showgraphics < 0)
        args->showgraphics = (args->max_iterations / 25) < 1 ? 1 : (args->max_iterations / 25);
    /* TODO
     * check for bounderies
     */


    return 0;
}

void static print_args(ui_args * args, struct argp argp) {
    /* TODO: Print input settings ONLY IF VERBOSE AND QUIET ALLOWS THAT */

    if (args->quiet)
        return;

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
    Angle of atack    : %.3f\n\
    Mach speed        : %.3f\n\
    Order             : %d\n\
    CFL               : %.3f\n\
    Max interations   : %d\n\
    Residual threshold: %.3e\n",
    args->angle, args->mach, args->order, args->cfl, args->max_iterations, args->nr_threashold);

    char tempd[128], tempg[128], tmp[64];
    sprintf(tmp, "%d ", args->showdetails);
    sprintf(tempd, "details each %siteration%s, ", args->showdetails > 1 ? tmp : "", args->showdetails > 1 ? "s" : "");
    sprintf(tmp, "%d ", args->showgraphics);
    sprintf(tempg, "graphics each %siteration%s, ", args->showgraphics > 1 ? tmp : "", args->showgraphics > 1 ? "s" : "");

    if (args->verbose || args->showdetails || args->showgraphics)
        printf("\n"BOLD"Verbose info: "RESET"%s%s%s%s\b\b \n",
    (/*!args->f_ds && !args->f_plot && !args->f_pc && !args->f_residue && */!args->verbose) ? "normal, " : "",
    args->verbose      ? "more verbose, " : "",
    args->showdetails  ? tempd : "",
    args->showgraphics ? tempg : ""/*,
    args->f_ds         ? "data structure, " : "",
    args->f_plot       ? "plot data, " : "",
    args->f_pc         ? "pressure coeficients (Cp), " : "",
    args->f_residue    ? "norm residue, " : ""*/);

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
    if (arg == NULL) {
        printf("\n%c -- invalid argument: '%s'\n\n", key, "");
        argp_usage(state);
    }
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
    if (arg == NULL) {
        printf("\n%c -- invalid argument: '%s'\n\n", key, "");
        argp_usage(state);
    }
    errno = 0;
    char * c;
    double val = strtod(arg, &c);
    if ((errno != 0) || (arg == c)) {
        printf("\n%c -- invalid argument: '%s'\n\n", key, arg);
        argp_usage(state);
    }

    return val;
}



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

/*
 *  This is just a set of macros to improve consistency in printing formats.
 *  It defines templates for:
 *      - INFO({printf}) : "INFO[FIL] FUNCTION: {printf} "
 *      - WARN({printf}) : 
 *
 *  This is also a pretty-printing bash colours. I mean, if enable, it will
 *  print using some colours.
 *
 */

/*                      TODO list

    * var size of MOD

    * * NOT FOR NOW / NOT IMPORTANT/RELEVANT **


    * * ALREADY DONE / PRETTY MUCH DONE **


*/



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                            CONFIGS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  INFO/WARN/ERR configs
 */


/*
 *  COLOURS configs
 */
#define USECOLOURS 1            /* if to use colours */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                          INFO/WARN/ERR
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define INFOMF(mod, fun, ...) \
                            printf("INFO[%3s] %s: ", mod, fun); \
                            printf(__VA_ARGS__); \
                            printf("\n")
#define INFO(...)           INFOMF(MOD, FUN, __VA_ARGS__)


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                            COLOURS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef USECOLOURS
    #define BOLD    "\033[1m"
    #define UND     "\033[4m"

    #define WARNT   "\033[1;33m"
    #define ERRT    "\033[1;31m"

    #define RESET   "\033[0m"
#else
    #define BOLD    ""
    #define UND     ""
    #define ERRT    ""
    #define ERRM    ""

    #define RESET   ""
#endif












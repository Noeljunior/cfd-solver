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
#define    _ICCS 3 /* the size of MOD */

/*
 *  COLOURS configs
 */
#define USECOLOURS 1            /* if to use colours */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                          INFO/WARN/ERR
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define MSGD(colour, type, mod, fun, ...)\
                         ({ char mi = mod != NULL && mod[0] != '\0' ? 1 : 0;\
                            char fi = fun != NULL && fun[0] != '\0' ? 1 : 0;\
                            printf(colour "%s"RESET"%s%*s%s %s%s", type,\
                                mi ? "[" : "", _ICCS, mi ? mod : "", mi ? "]" : "",\
                                fi ? fun : "", fi ? ": " : ""); \
                            printf(__VA_ARGS__); \
                            printf("\n");})

/* INFO messages */
#define INFOMF(...)         MSGD(INFOT, "INFO", MOD, FUN, __VA_ARGS__)
#define INFOM(...)          MSGD(INFOT, "INFO", MOD, "" , __VA_ARGS__)
#define INFO(...)           MSGD(INFOT, "INFO", "" , "" , __VA_ARGS__)

/* WARNING messages */
#define WARNMF(...)         MSGD(WARNT, "WARN", MOD, FUN, __VA_ARGS__)
#define WARNM(...)          MSGD(WARNT, "WARN", MOD, "" , __VA_ARGS__)
#define WARN(...)           MSGD(WARNT, "WARN", "" , "" , __VA_ARGS__)

/* ERROR messages */
#define ERRMF(...)          MSGD(ERRT,  "ERR" , MOD, FUN, __VA_ARGS__)
#define ERRM(...)           MSGD(ERRT,  "ERR" , MOD, "" , __VA_ARGS__)
#define ERR(...)            MSGD(ERRT,  "ERR" , "" , "" , __VA_ARGS__)

/* ERROR messages and exit */
#define AERRMF(...)         MSGD(ERRT,  "ERR" , MOD, FUN, __VA_ARGS__); exit(0);
#define AERRM(...)          MSGD(ERRT,  "ERR" , MOD, "" , __VA_ARGS__); exit(0);
#define AERR(...)           MSGD(ERRT,  "ERR" , "" , "" , __VA_ARGS__); exit(0);



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                            COLOURS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef USECOLOURS
    #define BOLD    "\033[1m"
    #define UND     "\033[4m"
    #define RED     "\033[31m"
    #define GREEN   "\033[32m"
    #define YELLOW  "\033[33m"
    #define BLUE    "\033[34m"
    #define PURPLE  "\033[35m"
    #define CYAN    "\033[36m"
    #define RESET   "\033[0m"


    #define INFOT   PURPLE BOLD

    #define WARNT   YELLOW BOLD

    #define ERRT    RED BOLD
#else
    #define BOLD    ""
    #define UND     ""
    #define RED     ""
    #define GREEN   ""
    #define YELLOW  ""
    #define BLUE    ""
    #define PURPLE  ""
    #define CYAN    ""
    #define ERRT    ""
    #define ERRM    ""

    #define RESET   ""
#endif













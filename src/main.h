#include "algorithm.h"

#define pi acos(-1.0)
#define PI 3.141592653589793238

typedef unsigned int uint;


/*
 * The global variables
 */
args *      inargs;             /* The input args */
mesh *      inmesh;             /* The input file */



/*
 * The funtions' headers
 */
args *      read_args(char **argv);
void        print_input(args * ina, mesh * inm);




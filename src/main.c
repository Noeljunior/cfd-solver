#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#include "main.h"


/*  TODO LIST:
    - add an option to select compatible mode or normal mode explicity
        . try to read the file/stdin and understand as default
    - add the ability to read configs as a file


*/

/* Print info related to this file */
#define INFO(i) printf("INFO[MA]: %s\n", i)

int main(int argc, char **argv) {
    #ifdef OLD                      /* The OLD args */
    assert(argc < 13, "Insufficient arguments given");
    INFO("Starting in compatible mode");
    #else                           /* The new args */
    assert(argc < 11, "Insufficient arguments given");
    INFO("Starting in normal mode");
    #endif
    
    inargs = read_args(argv);
    
    FILE * f = fopen(argv[1], "r");
    inmesh = ds_read_file(f);
    
    /*print_input(inargs, inmesh);/**/
    
    compute(inargs, inmesh);
    
    return 0;
}

/*
 * Read the input arguments
 */
args * read_args(char **argv) {
    args * in = (args *) malloc(sizeof(args));
    #ifdef OLD                      /* The OLD args */
    in->ORDER       = atoi(argv[ 2]);
    in->TP          = atoi(argv[ 4]);
    in->CFL         = atof(argv[ 5]);
    in->F_DS        = atoi(argv[ 6]);
    in->F_MP        = atoi(argv[ 8]);
    in->F_CP        = atoi(argv[ 9]);
    in->F_NR        = atoi(argv[10]);
    in->ITER_MAX    = atoi(argv[11]);
    in->NRT         = atof(argv[12]);
    #else                           /* The new args */
    in->ORDER       = atoi(argv[ 2]);
    in->CFL         = atof(argv[ 3]);
    in->ITER_MAX    = atoi(argv[ 4]);
    in->NRT         = atof(argv[ 5]);
    in->F_DS        = atoi(argv[ 6]);
    in->F_MP        = atoi(argv[ 7]);
    in->F_CP        = atoi(argv[ 8]);
    in->F_NR        = atoi(argv[ 9]);
    in->TP          = atoi(argv[10]);
    #endif
    return in;
}


void print_input(args * ina, mesh * inm) {
    /* printing the arguments */
    printf("%-26s: %d\n", "Order",                      ina->ORDER);
    printf("%-26s: %f\n", "CFL Condiction",             ina->CFL);
    printf("%-26s: %d\n", "Maximum iterations",         ina->ITER_MAX);
    printf("%-26s: %g\n", "Norm residue threshold",     ina->NRT);
    printf("%-26s: %s\n", "Print data structure",       ina->F_DS > 0 ? "yes" : "no" );
    printf("%-26s: %s\n", "Print data to plot",         ina->F_MP > 0 ? "yes" : "no" );
    printf("%-26s: %s\n", "Print pressure coeficients", ina->F_CP > 0 ? "yes" : "no" );
    printf("%-26s: %s\n", "Print norm residue",         ina->F_NR > 0 ? "yes" : "no" );
    printf("%-26s: %d\n", "Test problem",               ina->TP);
    
    /* Print the mesh */
    int i, c;
    
    printf("\nVertices (%d):\n", inm->novertices);
    for (i = 0, c = 0; i < inm->novertices; i++) {
        printf("[%10.4f, %10.4f]", inm->vertices[i].x, inm->vertices[i].y);
        if (++c >= 3) printf("\n", c = 0);
        else          printf("   ");
    }
    printf("\n");/**/
    
    printf("\nEdges (%d):\n", inm->noedges);
    for (i = 0, c = 0; i < inm->noedges; i++) {
        printf("[%6d, %6d; %3s]", inm->edges[i].a, inm->edges[i].b,
            inm->edges[i].border == INNER ? "IN" : "OUT");
        if (++c >= 3) printf("\n", c = 0);
        else          printf("        ");
    }
    printf("\n");
    
    printf("\nTriangles (%d):\n", inm->notriangles);
    for (i = 0, c = 0; i < inm->notriangles; i++) {
        printf("[%6d, %6d, %6d]", inm->triangles[i].a, inm->triangles[i].b, inm->triangles[i].c);
        if (++c >= 3) printf("\n", c = 0);
        else          printf("   ");
    }
    printf("\n");/**/
}













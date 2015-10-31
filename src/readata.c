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

#include "readata.h"

#include "infocolour.h"

#include <stdio.h>
#include <stdlib.h>

/*
 *  This file name, used in INFO/WARN/ERR msgs
 */
/*static char *MOD = "RD";*/


/*                      TODO list

    * detect file type/version and decide
    * option to read from stdin
    * a function to free this arrays
    * pointerify the struct
    * create a function to force read like specific format and source
        and a enum for type

*/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             FUNCTION DECLARATION and globals definitions
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             EXTERNAL INTERFACE
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
cfdrd_ds cfdrd_readfile_auto(char *path) {
    /*
     * This computes the CFD problem as Juan does
     */
    /*static char *FUN = "readfile_auto()";*/
    cfdrd_ds ds;


    double **vertices;
    int **edges, **triangles;
    int    sizev, sizee, sizet;

    /* JUAN file type */
    int i, t1, t2, t3;
    FILE * f = fopen(path, "r");
    char trash[1024];
    int  trashi;
    fscanf(f, "%s %d", trash, &trashi);           /* trash */
    fscanf(f, "%s %d", trash, &trashi);           /* trash */
    fscanf(f, "%s %d", trash, &sizev);    /* number of vertices */

    vertices       = (double **) malloc(sizeof(double *) * sizev);
    vertices[0]    = (double * ) calloc(sizev * 2, sizeof(double));
    for (i = 1; i < sizev; i++)
        vertices[i] = vertices[0] + i * 2;

    for (i = 0; i < sizev; i++) {           /* vertices */
        fscanf(f, "%lf %lf %d %d", &vertices[i][0], &vertices[i][1], &trashi, &trashi);
    }

    fscanf(f, "%s %d", trash, &sizee);       /* number of edges */

    edges       = (int **) malloc(sizeof(int *) * sizee);
    edges[0]    = (int * ) calloc(sizee * 3, sizeof(int));
    for (i = 1; i < sizee; i++)
        edges[i] = edges[0] + i * 3;

    for (i = 0; i < sizee; i++) {              /* edges */
        fscanf(f, "%d %d %d", &edges[i][0], &edges[i][1], &edges[i][2]);
        edges[i][0]--;
        edges[i][1]--;
        edges[i][2]--;
    }

    fscanf(f, "%s %d", trash, &sizet);   /* number of triangles */
    triangles       = (int **) malloc(sizeof(int *) * sizet);
    triangles[0]    = (int * ) calloc(sizet * 3, sizeof(int));
    for (i = 1; i < sizet; i++)
        triangles[i] = triangles[0] + i * 3;

    for (i = 0; i < sizet; i++) {         /* triangles */
        fscanf(f, "%d %d %d %d", &t1, &t2, &t3, &trashi);
        triangles[i][0]  = t1 - 1;
        triangles[i][1]  = t2 - 1;
        triangles[i][2]  = t3 - 1;
    }

    fscanf(f, "%s", trash);                       /* EOF */
    fclose(f);

    /* TODO read directly to here */
    ds.vertices  = vertices;
    ds.sizev     = sizev;

    ds.edges     = edges;
    ds.sizee     = sizee;

    ds.triangles = triangles;
    ds.sizet     = sizet;

    return ds;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             READATA
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


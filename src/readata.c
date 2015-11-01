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
static char *MOD = "RD";


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
void        read_new(FILE * f, char * fl, cfdrd_ds * ds);
void        read_old(FILE * f, char * fl, cfdrd_ds * ds);





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             EXTERNAL INTERFACE
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
cfdrd_ds * cfdrd_readfile_auto(char *path, int rstdin, int force, int bequiet) {
    static char *FUN = "cfdrd_readfile_auto()";
    cfdrd_ds * ds = (cfdrd_ds *) malloc(sizeof(cfdrd_ds));
    FILE *f;
    char *fl;
    size_t s;

    /* select file source */
    if (rstdin) f = stdin;
    else        f = fopen(path, "r");

    /* read the first line*/
    fl = NULL;
    getline(&fl, &s, f);

    /* call the righ parser */
    if      (force == 1) read_new(f, fl, ds);
    else if (force == 2) read_old(f, fl, ds);
    else {
        /* take a choice based on first line */
        int i, old = 1, new = 1;
        char *cmp = "MeshVersionFormatted";
        char tmp[1024];
        sscanf(fl, "%s", tmp);
        for (i = 0; i < 20; i++) {
            if (tmp[i] != cmp[i]) {
                old = 0;
                break;
            }
        }

        if (sscanf(fl, "%d", &i) != 1)
            new = 0;

        if (old == new) { /* oops! */
            AERRMF("Hey! Give me a well known mesh file, please...");
        }


        if      (new == 1) read_new(f, fl, ds);
        else if (old == 1) read_old(f, fl, ds);

        INFOMF("Input file was detected as %s file format", new ? "new" : "old");
    }
    free(fl);

    /* close the file if needed */
    if (!rstdin)
        fclose(f);

    return ds;
}

void cfdrd_free(cfdrd_ds * ds) {
    if (!ds->vertices) {
        free(ds->vertices[0]);
        free(ds->vertices);
    }
    if (!ds->edges) {
        free(ds->edges[0]);
        free(ds->edges);
    }
    if (!ds->triangles) {
        free(ds->triangles[0]);
        free(ds->triangles);
    }
    free(ds);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             READATA
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void read_new(FILE * f, char * fl, cfdrd_ds * ds) {
    /*
     *  standard file type
     */
    int i;

    sscanf(fl, "%d", &ds->sizev);              /* number of vertices */
    ds->vertices       = (double **) malloc(sizeof(double *) * ds->sizev);
    ds->vertices[0]    = (double * ) calloc(ds->sizev * 2, sizeof(double));
    for (i = 0; i < ds->sizev; i++) {           /* vertices */
        ds->vertices[i] = ds->vertices[0] + i * 2;
        fscanf(f, "%lf %lf", ds->vertices[i], ds->vertices[i] + 1);
    }

    fscanf(f, "%d", &ds->sizee);                 /* number of edges */
    ds->edges       = (int **) malloc(sizeof(int *) * ds->sizee);
    ds->edges[0]    = (int * ) calloc(ds->sizee * 3, sizeof(int));
    for (i = 0; i < ds->sizee; i++) {             /* edges */
        ds->edges[i] = ds->edges[0] + i * 3;
        fscanf(f, "%d %d %d", ds->edges[i], ds->edges[i] + 1, ds->edges[i] + 2);
    }

    fscanf(f, "%d", &ds->sizet);             /* number of triangles */
    ds->triangles       = (int **) malloc(sizeof(int *) * ds->sizet);
    ds->triangles[0]    = (int * ) calloc(ds->sizet * 3, sizeof(int));
    for (i = 0; i < ds->sizet; i++) {         /* triangles */
        ds->triangles[i] = ds->triangles[0] + i * 3;
        fscanf(f, "%d %d %d", ds->triangles[i], ds->triangles[i] + 1, ds->triangles[i] + 2);
    }
}

void read_old(FILE * f, char * fl, cfdrd_ds * ds) {
    /*
     *  JUAN file type
     */
    int i;
    //FILE * f = fopen(path, "r");
    char trash[1024];
    int  trashi;
    sscanf(fl, "%s %d", trash, &trashi);          /* trash */
    fscanf(f, "%s %d", trash, &trashi);           /* trash */

    fscanf(f, "%s %d", trash, &ds->sizev);    /* number of vertices */
    ds->vertices       = (double **) malloc(sizeof(double *) * ds->sizev);
    ds->vertices[0]    = (double * ) calloc(ds->sizev * 2, sizeof(double));
    for (i = 0; i < ds->sizev; i++) {           /* vertices */
        ds->vertices[i] = ds->vertices[0] + i * 2;
        fscanf(f, "%lf %lf %d %d", ds->vertices[i], ds->vertices[i] + 1, &trashi, &trashi);
    }

    fscanf(f, "%s %d", trash, &ds->sizee);       /* number of edges */
    ds->edges       = (int **) malloc(sizeof(int *) * ds->sizee);
    ds->edges[0]    = (int * ) calloc(ds->sizee * 3, sizeof(int));
    for (i = 0; i < ds->sizee; i++) {              /* edges */
        ds->edges[i] = ds->edges[0] + i * 3;
        fscanf(f, "%d %d %d", ds->edges[i], ds->edges[i] + 1, ds->edges[i] + 2);
        ds->edges[i][0]--;
        ds->edges[i][1]--;
        ds->edges[i][2]--;
    }

    fscanf(f, "%s %d", trash, &ds->sizet);   /* number of triangles */
    ds->triangles       = (int **) malloc(sizeof(int *) * ds->sizet);
    ds->triangles[0]    = (int * ) calloc(ds->sizet * 3, sizeof(int));
    for (i = 0; i < ds->sizet; i++) {         /* triangles */
        ds->triangles[i] = ds->triangles[0] + i * 3;
        fscanf(f, "%d %d %d %d", ds->triangles[i], ds->triangles[i] + 1, ds->triangles[i] + 2, &trashi);
        ds->triangles[i][0]--;
        ds->triangles[i][1]--;
        ds->triangles[i][2]--;
    }

    fscanf(f, "%s", trash);                       /* EOF */
}








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
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>

/*
 *  This file name, used in INFO/WARN/ERR msgs
 */
static char *MOD = "RDP";
static char *COL = CYAN;


#define MAXPFUNCS   16
#define MAXPARAMS   16


/*                      TODO list

    * ability to convert from old to new file format


    * * NOT FOR NOW / NOT IMPORTANT/RELEVANT **


    * * ALREADY DONE / PRETTY MUCH DONE **
    * a function to free this arrays
    * pointerify the struct
    * create a function to force read like specific format and source
        and a enum for type
    * detect file type/version and decide
    * option to read from stdin

*/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             FUNCTION DECLARATION and globals definitions
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void        parse_cfdm_file(FILE * f, char * fl, cfdrd_ds * ds, int bequiet);
void        parse_old_file(FILE * f, char * fl, cfdrd_ds * ds, int bequiet);
double *    radius_spline(int N ,double *tt ,double *xx, double *yy);



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             EXTERNAL INTERFACE
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
cfdrd_ds * cfdrd_readfile_auto(char *path, int rstdin, int force, int bequiet) {
    static char *FUN = "cfdrd_readfile_auto()";
    cfdrd_ds * ds = (cfdrd_ds *) calloc(1, sizeof(cfdrd_ds));
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
    if      (force == 1) parse_cfdm_file(f, fl, ds, bequiet);
    else if (force == 2) parse_old_file(f, fl, ds, bequiet);
    else {
        /* take a choice based on first line */
        int i, old = 1, cfdm = 1;
        char *cmp = "MeshVersionFormatted";
        char tmp[1024];
        sscanf(fl, "%s", tmp);
        for (i = 0; i < 20; i++) {
            if (tmp[i] != cmp[i]) {
                old = 0;
                break;
            }
        }

        if (old == cfdm) { /* oops! */
            AERRMF("Hey! Give me a well known mesh file, please...");
        }

        if (!bequiet)
            INFOMF("Input file was detected as %s file format", cfdm ? "cfdm" : "old");

        if      (cfdm == 1) parse_cfdm_file(f, fl, ds, bequiet);
        else if (old == 1) parse_old_file(f, fl, ds, bequiet);
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
    /*if (!ds->edges) {
        free(ds->edges[0]);
        free(ds->edges);
    }*/

    if (!ds->triangles) {
        free(ds->triangles[0]);
        free(ds->triangles);
    }

    if (ds->sizep && !ds->params) {
        int i;
        for (i = 0; i < ds->sizep; i++)
            free(ds->params[i]);
        free(ds->params);
    }

    free(ds);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             READATA
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void parse_cfdm_file(FILE * f, char * fl, cfdrd_ds * ds, int bequiet) {
    static char *FUN = "cfdrd_parse_input_file()";
    //cfdrd_ds * ds = (cfdrd_ds *) calloc(1, sizeof(cfdrd_ds));

    //FILE   *f;
    size_t  s = strlen(fl);
    char   *l = (char *) malloc(sizeof(char) * s);
    strncpy(l, fl, s);

    /*if (rstdin) f = stdin;
    else        f = fopen(path, "r");*/

    ssize_t       rc = s;
    int           state = 0;
    int           rv    = 0;
    unsigned long cl    = 0;
    int           cv    = 0,
                  ct    = 0;

    int i;

    do {
        ++cl;
        /* empty lines are just ignored */
        if (l[0] == '\n')   continue;

        /* comment line are just ignored */
        if (l[0] == '#')   continue;

        l[rc-1] = '\0';
        /* a to-print line */
        if (l[0] == '$') {
            if (!bequiet)
                INFOMF("at line %lu: %s", cl, l + 1);
            continue;
        }

        if (state == 0) {
            /* at this state we expect to read the line containing sizes, otherwise abort */
            ds->chord = -1;
            rv = sscanf(l, "%d %d %lf", &ds->sizev, &ds->sizet, &ds->chord);
            if (rv < 2)
                AERRMF("at line %lu: couldn't understand how many vertices and triangles are - '%s'", cl, l);

            /* alloc mem */
            ds->vertices       = (double **) malloc(sizeof(double *) * ds->sizev);
            ds->vertices[0]    = (double * ) calloc(ds->sizev * 4, sizeof(double));
            for (i = 0; i < ds->sizev; i++)  ds->vertices[i] = ds->vertices[0] + i * 4;

            /*ds->border = (int *)    calloc(ds->sizev, sizeof(int));
            ds->radius = (double *) calloc(ds->sizev, sizeof(double));*/

            ds->triangles       = (int **) malloc(sizeof(int *) * ds->sizet);
            ds->triangles[0]    = (int * ) calloc(ds->sizet * 3, sizeof(int));
            for (i = 0; i < ds->sizet; i++) ds->triangles[i] = ds->triangles[0] + i * 3;

            state = 1;
        } else if (state > 0) {
            /* by now we should be accepting line starting with 'params'
               if a line do not start with that, we assume that are either vertices or triangles lines */

            if (!strncmp("params ", l, rc < 7 ? rc : 7)) {
                int border = -1, func = -1, args = -1, len = -1;
                rv = sscanf(l+7, "%d %d %d%n", &border, &func, &args, &len);

                if (border < 0 || func < 0 || args < 1)
                    AERRMF("at line %lu: something is wrong with these parameters - '%s'", cl, l);

                ds->sizep++;
                ds->params = (double **) realloc(ds->params, sizeof(double *) * ds->sizep);
                ds->params[ds->sizep-1] = (double *) malloc(sizeof(double) * (3 + args));
                ds->params[ds->sizep-1][0] = border;
                ds->params[ds->sizep-1][1] = func;
                ds->params[ds->sizep-1][2] = args;

                int n = 0;
                for (i = 0; i < args; i++) {
                    n = 0;
                    rv = sscanf(l+7+len+n, "%lf%n", ds->params[ds->sizep-1] + 3 + i, &n);
                    len += n;
                    if (rv != 1)
                        AERRMF("at line %lu: something is wrong with these parameters - '%s'", cl, l);
                }

                INFOMF("at line %lu: new params read", cl);
            } else if (state == 1) { /* reading a vertice */
                double x, y, b, r;
                b = -1.0; r = 0.0;
                rv = sscanf(l, "%lf %lf %lf %lf", &x, &y, &b, &r);
                if (rv < 2)
                    AERRMF("at line %lu: something is wrong with this vertice - '%s'", cl, l);

                ds->vertices[cv][0] = x;
                ds->vertices[cv][1] = y;
                ds->vertices[cv][2] = b;
                ds->vertices[cv][3] = r;
                /*ds->border[cv] = b;
                ds->radius[cv] = r;*/

                if (++cv == ds->sizev)
                    state = 2;

                /* check : detect a trailling point :  */
                /*if (!isfinite(r))
                    printf("foud the trailing ;)\n");*/
            } else if (state == 2) { /* reading a triangle */
                int v1, v2, v3;
                rv = sscanf(l, "%d %d %d", &v1, &v2, &v3);
                if (rv < 3)
                    AERRMF("at line %lu: something is wrong with this triangle - '%s'", cl, l);

                ds->triangles[ct][0] = v1;
                ds->triangles[ct][1] = v2;
                ds->triangles[ct][2] = v3;

                if (++ct == ds->sizet)
                    state = -1;
            } else break;
        }

    } while (((rc = getline(&l, &s, f)) != -1) && state >= 0);

    free(l);
}

void parse_old_file(FILE * f, char * fl, cfdrd_ds * ds, int bequiet) {
    /*
     *  JUAN file type
     */
    int i;

    double      **vertices;
    int         sizev;
    int         **edges;
    int         sizee;
    int         **triangles;
    int         sizet;

    char trash[1024];
    int  trashi;
    sscanf(fl, "%s %d", trash, &trashi);          /* trash */
    fscanf(f, "%s %d", trash, &trashi);           /* trash */

    fscanf(f, "%s %d", trash, &sizev);    /* number of vertices */
    vertices       = (double **) malloc(sizeof(double *) * sizev);
    vertices[0]    = (double * ) calloc(sizev * 4, sizeof(double));
    for (i = 0; i < sizev; i++) {           /* vertices */
        vertices[i] = vertices[0] + i * 4;
        fscanf(f, "%lf %lf %d %d", vertices[i], vertices[i] + 1, &trashi, &trashi);
        vertices[i][2] = -1;
    }

    fscanf(f, "%s %d", trash, &sizee);       /* number of edges */
    edges       = (int **) malloc(sizeof(int *) * sizee);
    edges[0]    = (int * ) calloc(sizee * 3, sizeof(int));
    for (i = 0; i < sizee; i++) {              /* edges */
        edges[i] = edges[0] + i * 3;
        fscanf(f, "%d %d %d", edges[i], edges[i] + 1, edges[i] + 2);
        edges[i][0]--;
        edges[i][1]--;
        edges[i][2]--;
    }

    fscanf(f, "%s %d", trash, &sizet);   /* number of triangles */
    triangles       = (int **) malloc(sizeof(int *) * sizet);
    triangles[0]    = (int * ) calloc(sizet * 3, sizeof(int));
    for (i = 0; i < sizet; i++) {         /* triangles */
        triangles[i] = triangles[0] + i * 3;
        fscanf(f, "%d %d %d %d", triangles[i], triangles[i] + 1, triangles[i] + 2, &trashi);
        triangles[i][0]--;
        triangles[i][1]--;
        triangles[i][2]--;
    }

    fscanf(f, "%s", trash);                       /* EOF */

    /* compute the borders' radius */

    int sizes[2] = {0}; /* hardcoded maximum borders */
    for (i = 0; i < sizee; i++) {
        sizes[edges[i][2]]++;
    }

    /* airfoil border */
    double *tt = malloc(sizeof(double) * (sizes[0] + 1)),
           *xx = malloc(sizeof(double) * (sizes[0] + 1)), /* xx values must be strictly increasing */
           *yy = malloc(sizeof(double) * (sizes[0] + 1)),
           *curve_radius;

    tt[0] = 0;
    xx[0] = vertices[edges[0][0]][0];
    yy[0] = vertices[edges[0][0]][1];
    for (i = 1; i < sizes[0] + 1; i++) {
        tt[i] = i;
        xx[i] = vertices[edges[i - 1][1]][0];
        yy[i] = vertices[edges[i - 1][1]][1];
    }

    curve_radius = radius_spline(sizes[0] + 1, tt, xx, yy);

    for (i = 0; i < sizes[0]; i++) {
        vertices[edges[i][0]][2] = 0;
        if (!isfinite(curve_radius[i]))
            continue;
        vertices[edges[i][0]][3] = curve_radius[i];
    }

    /* far-end border */
    printf("\n");
    double sum = 0.0;
    for (i = 0; i < sizes[1]; i++) {
        sum += sqrt(pow(vertices[edges[i + sizes[0]][0]][0], 2) + pow(vertices[edges[i + sizes[0]][0]][1], 2));
    }
    sum /= sizes[1];
    for (i = 0; i < sizes[1]; i++) {
        vertices[edges[i + sizes[0]][0]][2] = 1;
        vertices[edges[i + sizes[0]][0]][3] = sum;
    }

    /* chord */
    /*double min = DBL_MAX, max = DBL_MIN;
    int imin, imax;
    for (i = 0; i < sizes[0]; i++) {
        if (vertices[edges[i][0]][0] < min) {
            min = vertices[edges[i][0]][0] ;
            imin = edges[i][0];
        }
        if (vertices[edges[i][0]][0] > max) {
            max = vertices[edges[i][0]][0] ;
            imax = edges[i][0];
        }
    }

    ds->chord = sqrt(pow(ds.vertices[imin][0] - ds.vertices[imax][0], 2) + pow(ds.vertices[imin][1] - ds.vertices[imax][1], 2));
    */ds->chord = 1.0;


    /* Free! */
    free(tt);
    free(xx);
    free(yy);
    free(curve_radius);
    free(edges[0]);
    free(edges);

    /* put this data into ds */
    ds->vertices  = vertices;
    ds->sizev     = sizev;
    ds->triangles = triangles;
    ds->sizet     = sizet;
}

double * radius_spline(int N ,double *tt ,double *xx, double *yy) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc () ;
    gsl_spline *spline_x = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline *spline_y = gsl_spline_alloc (gsl_interp_cspline, N);

    double *curve_radius = malloc(sizeof(double) * N);

    gsl_spline_init(spline_x, tt, xx, N);
    gsl_spline_init(spline_y, tt, yy, N);

    int i ; 
    double dxi, dyi, d2xi, d2yi;

    for (i = 0 ; i < N; i++) {
        dxi  = gsl_spline_eval_deriv( spline_x, tt[i], acc);
        dyi  = gsl_spline_eval_deriv( spline_y, tt[i], acc);
        d2xi = gsl_spline_eval_deriv2(spline_x, tt[i], acc);
        d2yi = gsl_spline_eval_deriv2(spline_y, tt[i], acc);
        curve_radius[i] = pow (dxi * dxi + dyi * dyi, 1.5) / (dxi * d2yi - dyi*d2xi);
    }

    gsl_spline_free(spline_x);
    gsl_spline_free(spline_y);
    gsl_interp_accel_free(acc);

    return curve_radius;
}




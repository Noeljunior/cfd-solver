#include <stdio.h>

/* Versions of input files */
#define VERSION_G1 1
#define VERSION_J1 2

typedef unsigned int uint;


/*
 * The input arguments
 */
#ifndef DS_IFDEF
#define DS_IFDEF
struct args {
    int     ORDER;              /* Reconstruction order */
    double  CFL;                /* CFL Condiction */
    int     ITER_MAX;           /* Maximum iterations */
    double  NRT;                /* Norm residue threshold */
    int     F_DS;               /* FLAG: print data structure */
    int     F_MP;               /* FLAG: print data to plot */
    int     F_CP;               /* FLAG: print pressure coeficients, Cp */
    int     F_NR;               /* FLAG: print norm residue */
    int     TP;                 /* Test problem */

    int     NVERTICES;          /* Number of vertices */

};

/*
 * The input file content
 */
typedef struct mesh {
    int     dimension;          /* Dimension */

    int     novertices;         /* Number of vertices */
    struct  coord *vertices;    /* The vertices */

    int     noedges;            /* Number of edges */
    struct  edge *edges;        /* The edges */

    int     notriangles;        /* Number of triangles */
    struct  triangle *triangles;/* The triangles */

    int     novt;               /* Number of vertices plus triangles */
    int     noc;                /* Number of coeficients of Taylor 2D polynomial */
    int     nommt;              /* Number of gauss points necessary to the MMT momemtum calcs */
    double ** mmt;              /* Momentus matrix */
    int     lpd;                /* Legendre's polynomial degree */
    int     Aexp;               /* Exponent of the weights of A matrix */
    int     norv;               /* Number of reconstruction variables */
    double  kdelta;             /* Auxiliary constant of limiters */
    double  m1, m2;             /* Mach's numbers used in limiters calculus */
    struct  face *faces;        /* The faces TODO understand wtf is this */

    double  *gps[3];            /* Gauss points in [-1, 1] */
    double  *gws[3];            /* Gauss weights */
    struct  coord *tcentroids;  /* Triangles' centroid coords */

} mesh;


/*** A 2D coordinate ***/
typedef struct coord {
    double  x, y;
} coord;

/*** An edge ***/
typedef enum border { INNER = 0, OUTTER = 1 } border;
typedef struct edge {
    uint    a, b;
    border  border;
} edge;

/*** A triangle ***/
typedef struct triangle {
    uint    a, b, c;
} triangle;

/*** A face ***/
typedef struct face {
    uint    vi, vf,
            tl, tr,
            vl, vr;
} face;

#endif

typedef struct args args;

/*
 * The funtions' headers
 */

mesh *      ds_read_file(FILE *fid);



void        abortc(char *why, int errno);
void        abortw(char *why);
void        aborte();
void        assert(int cond, char *why);
void        nassert(int cond, char *why);

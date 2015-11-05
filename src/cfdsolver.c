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

#include "cfdsolver.h"

#include "infocolour.h"

#include "meshviewer.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*
 *  This file name, used in INFO/WARN/ERR msgs
 */
static char *MOD = "SLV";
static char *COL = BLUE;

/*                      TODO list
    * inflow angle and mach_inflow as input args
    * -d[detailed output of simulation] parameter: ability to chose the frequency of output
    * -g[raphics]   show meshviewer at a given frequency



    * receive the parametric naca function and cl and cd references


    * IMPLEMENT the stencil!

    * Print drag and lift coeficients


    * solver function RETURNS SOME CLASSIFICATION

    * ability to save stats to a file!
        * choose what to output

    * build a output visualizer

    * integrate this with the meshviewer!
        + pressures?
    * only show the mesh if asked by user by argument
    * stop signal handling after solver finish

    * REFACTOR MESHVIEWER



    * * NOT FOR NOW / NOT IMPORTANT/RELEVANT **
    + verbose to files for stats/plots
    * configurable number of digits of printed residue/objective_value
    * option mesh draw



    * * ALREADY DONE / PRETTY MUCH DONE **
    * fix mutexes deadlock
    + INFO
        INFO + WARN + ERR based on verbosity
        bold/colours CONFIGURABLE
    * ifdef for debugs
    +   verbosity/quiet
        quiet < normal/default < verbose
        force output the objective value
    * check what args should be given
    * TIMER align time prints by max size of seconds (like number of iter)
    * put the INFO/WARN/ERR in a separated .h
    * add ability to count time of init functions
    * make the edge reading a like a face building
    * percorrer todas as faces fronteira de forma seemless
    * Percorrer todos os vertices fronteira de forma seemless

*/


/*
 *  Marcros Helpers
 */
#define forvar(var, s, cond) for (var = s; var < (cond); var++)

#define fori(cond)      forvar(i, 0, cond)
#define forj(cond)      forvar(j, 0, cond)
#define fork(cond)      forvar(k, 0, cond)
#define forn(cond)      forvar(n, 0, cond)
#define form(cond)      forvar(m, 0, cond)

#define forin(n, cond)  forvar(i, n, cond)
#define forjn(n, cond)  forvar(j, n, cond)
#define forkn(n, cond)  forvar(k, n, cond)
#define fornn(m, cond)  forvar(n, m, cond)
#define formn(n, cond)  forvar(m, n, cond)

#define mdouble(n)      ((double *)  malloc(sizeof(double) * (n)))
#define mdoublep(n)     ((double **) malloc(sizeof(double*) * (n)))
#define mint(n)         ((int *)     malloc(sizeof(int) * (n)))

#define cdouble(n)      ((double *)  calloc(n, sizeof(double)))
#define cint(n)         ((int *)     calloc(n, sizeof(int)))

#define pow2(b)         ((b)*(b))
#define pow3(b)         ((b)*(b)*(b))
#define l2dist(a, b)    sqrt(pow2((a)) + pow2((b)))

#define zs(var, arr, s) for(var = 0; var < s; var++) arr[var] = 0

#define indexof(a, i)   ((int) ((i) - (a)))

#define PI              3.14159265358979323846

/* how many decimal cases to print */
#define ITRESSIZE       5
#define FRESSIZE        15


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             DATA STRUCTURE and typedefs
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define NOU 4                   /* Number of reconstruction vars \Tu/ */

typedef struct cfds_mesh mesh;

typedef unsigned int uint;

typedef enum border { NONE = -1, INNER = 0, OUTTER = 1 } border;
typedef enum flow { INTERNAL = 0, WALL = 1, INFLOW = 2, OUTFLOW = 3 } flow;
typedef enum curved { NN = 0, AB = 1, BC = 2, CA = 3 } curved;

/*
 * The input file content
 */
struct cfds_mesh {
                                /* vertices, triangles and faces */
    int             novertices;             /* Number of vertices */
    struct          vertex *vertices;       /* The vertices */

    int             notriangles;            /* Number of triangles */
    struct          triangle *triangles;    /* The triangles */

    int             nofaces;                /* Number of faces */
    struct          face *faces;            /* The faces */

    int             noborders;              /* number of different borders */
    int             *noedgeface;            /* no of faces in each boder */
    struct          face **edgeface;        /* a pointer to the first element of each border */

                                /* problem settings */
    int             order;                  /* Reconstruction order */
    int             nommt;                  /* Number of gauss points necessary to the MMT momemtum calcs */
    int             notc;                   /* Number of coeficients of Taylor 2D polynomial \NC/ */
    int             lpd;                    /* Legendre's polynomial degree \GPL/ */
    double          cfl;                    /* CFL condiction \CFL/ */
    int             max_iter;               /* Max iterations to stop simulation \ITERMAX/ */
    double          nrt;                    /* Threshold to stop simulation \TOL/ */

                                /* constants */
    double          inflow_angle,           /* Inflow angle */
                    gamma,                  /* air constant */
                    r,                      /* ideal gas costant */
                    mach_inflow,            /* mach number in inflow */
                    t_inflow,               /* inflow's temperature */
                    u_inflow[NOU],          /* inflow primitives */
                    tt_inflow,              /* total temperature */
                    pt_inflow,              /* total pressure */
                    q_inflow,               /* |velocity| */
                    chord,                  /* airfoil chord */
                    cl_ref,                 /* reference lift coeficient */
                    cd_ref,                 /* reference drag coeficient */
                    k_delta;                /* Limiters auxiliary constant */
    int             Aexp;                   /* Exponent of the weights of A matrix \POT/ */
    double          M1, M2;                 /* Mach numbers @ limiters */

                                /* helpers */
    int             max_cvno;               /* the size of the biggest stencil */
    double          *gauss_points[3];       /* Gauss quadrature points */
    double          *gauss_weights[3];      /* Gauss quadrature weights */
    double          *gps[3];                /* Gauss points in [-1, 1] */
    double          *gws[3];                /* Gauss weights */

                                /* simulation */
    double          **uconserv,             /*  */
                    **u,                    /*  */
                    ***coef,                /*  */
                    **phi_chapel,           /*  */
                    **res;                  /*  */


                                /* settings */
    char            verbosity;              /* verbosity; 0: norma; 1: more verbose */
    char            quiet;                  /* be quiet all the time */
    int             showdetails;            /* if to show results each <showdetails> rungekutta iteration */
    int             showgraphics;           /* if to show graphics each <showgraphics> rungekutta iteration */
    char            fclassify;              /* if quiet is set, show the classification */
    char            plotviewer;             /* if to plot the meshviewer */

                                /* meshviewer */
    int             mesh_plot;              /* the identifier of plot in meshviewer */
    int             pressure_plot;          /* the identifier of plot in meshviewer */
    int             pressure_plotlog;       /* the identifier of plot in meshviewer (log scale) */
    int             velocity_plot;          /* the identifier of the velocity plot in meshviewer */
};


/*** A 2D coordinate ***/
typedef struct coord {
    double          x, y;
} coord;

/*** A 2D vertex ***/
typedef struct vertex {
    double          x, y;                   /* vertex' coordiantes */
    border          border;                 /* vertex borter type */
    flow            flow;                   /* vertex flow condiction type */
    double          curve_radius;           /* the curve radius, if in a border */
    double          *momentum;              /* the momemtum */
    coord           projection;             /* the projection */
    struct face     *preface, *posface;     /* when this vertex is in a border, the two adjacent faces also in border */

                            /* AS control volume */
    double          cva;                    /* Control volume area */
    int             cv_no;                  /* number of adjacent faces \NEE/ */
    //int             cv_no_1l;               /* TODO deprecated since we do not construct a stencil \NE1C/ */
    double          *m;                     /* gauss eliminations helper (?) \M[]/ */
    double          **pinv;                 /* Pseudoinverse \Pseudoinv_ULSP/ */

    double          a_ghostp;               /* Auxiliary var for ghost pressure calc \Aux_PressgVP/ */
    double          a_delta[NOU];           /* Auxiliary var \Aux_Delta/ */
    double          ghostp[NOU];            /* ghost pressure calc (limiters) \ugVP/ */
    struct vertex   **cv_stencil;

    int             aino;                   /* number of adjacent NON-BORDER faces \NumeroFIV/ */
    struct face     **aiface;               /* adjacent NON-BORDER faces \FIV/ */
} vertex;

/*** A triangle ***/
typedef struct triangle {
    vertex          *a, *b, *c;             /* the vertices that make this triangle */
    curved          curved;                 /* which pair of vertices are in border and curved, if any */
    struct face     *curvedface;            /* if curved, a pointer to that curved face */
    coord           centroid;               /* the centroid */
    border          border;                 /* in which border is this triangle */
} triangle;

/*** A face ***/
typedef struct face {
    vertex          *vi, *vf,               /* initial and final vertices \FACE:1:2/ */
                    *vl, *vr;               /* left and right vertices    \FACE:5:6/ */
    triangle        *tl, *tr;               /* left and right triangles   \FACE:3:4/ */
    double          curve_radius;           /* if this is curved, the curvature radius */
    coord           middle_point;           /* the middle point used when a face is in an edge */
    border          border;                 /* in which border is this face */

    coord           *gauss_points;          /* gauss points of this face  \XG/ */
    double          *gauss_weights;         /* gauss weights of this face \WG/ */
    coord           *normals;               /* normals of this face       \NXG/ */
    flow            *flow;                  /* type of flow flowing through each gausses' \tipoPG/ */
} face;

/*
 * Time measurement
 */
#define TIMEDESC 1024
typedef struct times {
    struct ttick    *tticks;                /* the time measurements */
    int             size;                   /* the size of above array */
    int             count;                  /* how many used untill now */
} times;

typedef struct ttick {
    char            desc[TIMEDESC];         /* a description of the this measurement */
    double          ticku;                  /* the actual time tick: user time */
    double          ticks;                  /* the actual time tick: sys  time */
    char            skip;                   /* -1: init; 0: tick; 1: skip the measurement */
} ttick;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             FUNCTION DECLARATION and globals definitions
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Solver fases
 */
void            compute(mesh * inm);
void            compute_faces(mesh * inm, int ** edges, int sizee);
void            compute_radius(mesh * inm);
void            compute_middle_points(mesh * inm);
void            compute_legendre(mesh * inm);
void            compute_centroids(mesh * inm);
void            compute_momentum(mesh * inm);
void            compute_projection(mesh * inm);
void            compute_stencil(mesh * inm);
void            compute_gausses(mesh * inm);
void            compute_coefmat_pseudoinverse(mesh * inm);
void            compute_prelimiters(mesh * inm);
void            compute_rungekutta5(mesh * inm);


/*
 * Solver helpers
 */
static void     compute_rk_convert(mesh * inm, double **uc, double **u);
static void     compute_rk_reconstruction(mesh * inm, double **u, double ***coef);
static double   compute_rk_polinomial(mesh * inm, double ***coef, int p, int k, double x, double y);
static void     compute_rk_limiters(mesh * inm, double **u, double ***coef, double **phi_chapel);
static double   compute_rk_polinomial_lim(mesh * inm, double **u , double ***coef, double **phi_chapel, int p, int k, double x, double y);
static void     compute_rk_r_primitive_lim(mesh * inm, double *vp, double **u , double ***coef, double **phi_chapel, int p, double x, double y);
static void     compute_rk_r_subinflow(mesh * inm, double * f, double press, double nx, double ny);
static void     compute_rk_r_suboutflow(mesh * inm, double * f, double * vp, double nx, double ny);
static void     compute_rk_r_roe(mesh * inm, double *F, double *vp1, double *vp2, double nx, double ny);
static void     compute_rk_residue(mesh * inm, double **u, double ***coef, double **phi_chapel, double **r, int iter);
void            compute_momentum_vertex(double * m, double ** gp, double ** gw, int o, int nommt, int t, coord vr, coord vi, coord vf, double r);
double *        compute_radius_spline(int N ,double *tt ,double *xx, double *yy);
flow            compute_wall_condiction(border b, double inflow_angle, double xx, double yy, double Nxx, double Nyy);

/*
 * meshviewer integration
 */
void            v_draw_rawmesh(mesh * inm);
void            v_draw_coefs(mesh * inm, int uselog);


/*
 * Time measurement
 */
times *         times_init(int maxmes);
void            times_tick(times * t, char * desc);
void            times_zero(times * t);
void            times_print(times * t, int verbose);
/* keep time-tracking */
times *         timemes;
int             supress_compute  = sizeof("compute");


/*
 * Signal handling
 */
void            sig_handler(int signo);
int             interrupt_solver = 0;
int             running_solver   = 0;

/*
 * Debug
 */
void            print2d(char * p, double **u, int l, int c, int new);
void            print3d(double ***u, int m, int l, int c, int new);



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             EXTERNAL INTERFACE
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
cfds_mesh * cfds_init(cfds_args * ina, double ** vertices, int sizev, int ** edges, int sizee, int ** triangles, int sizet) {
    /*
     * This initializes the internal data struct and returns it as ANONYMOUS
     */
    static char *FUN = "init()";

    /* check for input bounderies */
    


    /* new mesh data struct */
    mesh * inm = (mesh *) calloc(1, sizeof(mesh));

    /* copy input arguments to the internal data struct */
    inm->order                  = ina->order;
    inm->cfl                    = ina->cfl;
    inm->max_iter               = ina->max_iterations;
    inm->nrt                    = ina->nr_threashold;
    inm->inflow_angle           = ina->angle * PI / 180.0;
    inm->mach_inflow            = ina->mach;

    /* compute the values and set the constants */
    inm->nommt                  = (int) (ceil((inm->order + 1.0) / 2.0));
    inm->notc                   = (inm->order * (inm->order + 1)) / 2;
    inm->lpd                    = (int) (ceil(inm->order / 2.0));

    inm->gamma                  = 1.4;
    inm->r                      = 287.05;
    inm->t_inflow               = 233.0;
    inm->u_inflow[3]            = 25000.0;
    inm->tt_inflow              = inm->t_inflow * (1.0 + (inm->gamma - 1.0) / 2.0 * inm->mach_inflow * inm->mach_inflow);
    inm->pt_inflow              = inm->u_inflow[3] * pow(1.0 + (inm->gamma - 1.0) / 2.0 * inm->mach_inflow * inm->mach_inflow, inm->gamma / (inm->gamma - 1.0));
    inm->q_inflow               = inm->mach_inflow * sqrt(inm->gamma * inm->r * inm->t_inflow);
    inm->u_inflow[0]            = inm->u_inflow[3] / (inm->r * inm->t_inflow);
    inm->u_inflow[1]            = inm->q_inflow * cos(inm->inflow_angle);
    inm->u_inflow[2]            = inm->q_inflow * sin(inm->inflow_angle);
    inm->chord                  = 1.0;
    inm->cl_ref                 = 0.35169;
    inm->cd_ref                 = 0.022628;
    inm->k_delta                = 0.25;

    inm->Aexp                   = 1;
    inm->M1                     = 0.8;
    inm->M2                     = 0.85;

    /* copy settings */
    inm->verbosity              = ina->verbose;
    inm->quiet                  = ina->quiet;
    inm->showdetails            = ina->showdetails;
    inm->showgraphics           = ina->showgraphics;
    inm->fclassify              = ina->fclassify;

    /* initialize plot identifiers */
    inm->mesh_plot              = -1;
    inm->pressure_plot          = -1;
    inm->pressure_plotlog       = -1;
    inm->velocity_plot          = -1;

    /* copy vertices, edges and triangles to the internal data struct */
    int i;
    timemes = times_init(16);

    inm->novertices = sizev;
    inm->vertices = (vertex *) calloc(inm->novertices, sizeof(vertex));
    fori (sizev) {
        inm->vertices[i].x = vertices[i][0];
        inm->vertices[i].y = vertices[i][1];
        inm->vertices[i].border = NONE;
    }

    inm->notriangles = sizet;
    inm->triangles = (triangle *) calloc(inm->notriangles, sizeof(triangle));
    fori (sizet) {
        inm->triangles[i].a  = inm->vertices + triangles[i][0];
        inm->triangles[i].b  = inm->vertices + triangles[i][1];
        inm->triangles[i].c  = inm->vertices + triangles[i][2];
        inm->triangles[i].border = NONE;
        inm->triangles[i].curved = NN;
    }

    if (!inm->quiet) INFOMF("%d vertices, %d triangles and %d edges successfully added.", sizev, sizet, sizee);
    times_tick(timemes, "cfds_init()");

    /* compute the faces */
    compute_faces(inm, edges, sizee);
    times_tick(timemes, "compute_faces()" + supress_compute);

    /* VISUALIZER */
    if (inm->showgraphics) {
        mv_start(3, inm->quiet ? -1 : inm->verbosity);
        v_draw_rawmesh(inm);
    }

    return inm;
}


void cfds_solve(cfds_mesh * inm) {
    /*
     * This computes the CFD problem as Juan does
     */
    static char *FUN = "solve()";

    /* handle SINGINT so we can interrupt the solver; TODO maybe using a USERSIG is a better idea! */
    signal(SIGINT, sig_handler);
    /*INFOMF("Couldn't handle SIGINT signal. Not sure why (are you handling it before?)");*/

    /* Start count time now */
    times_zero(timemes);

    if (!inm->quiet && inm->verbosity) INFOMF("will cook every single var out there");

    compute_radius(inm);
    times_tick(timemes, "compute_radius()" + supress_compute);

    compute_middle_points(inm);
    times_tick(timemes, "compute_middle_points()" + supress_compute);

    compute_legendre(inm);
    times_tick(timemes, "compute_legendre()" + supress_compute);

    compute_centroids(inm);
    times_tick(timemes, "compute_centroids()" + supress_compute);

    compute_momentum(inm);
    times_tick(timemes, "compute_momentum()" + supress_compute);

    compute_projection(inm);
    times_tick(timemes, "compute_projection()" + supress_compute);

    compute_stencil(inm);
    times_tick(timemes, "compute_stencil()" + supress_compute);

    compute_gausses(inm);
    times_tick(timemes, "compute_gausses()" + supress_compute);

    compute_coefmat_pseudoinverse(inm);
    times_tick(timemes, "compute_coefmat_pseudoinverse()" + supress_compute);

    compute_prelimiters(inm);
    times_tick(timemes, "compute_prelimiters()" + supress_compute);

    if (!inm->quiet && inm->verbosity) INFOMF("everything is cooked... will now eat the CFD problem (runge kutta: 5 stages)");

    compute_rungekutta5(inm);
    times_tick(timemes, "compute_rungekutta5()" + supress_compute);

    /* Print time measurements */
    if (!inm->quiet) {
        if (inm->verbosity) printf("\n");
        times_print(timemes, inm->verbosity);
    }

    /* VISUALIZER */
    mv_wait();
}

void cfds_free(cfds_mesh * inm) {
    /*
     * This frees up every single bit allocated
     */
    static char *FUN = "free()";

    if (!inm->quiet && inm->verbosity)
        WARNMF("[NIY-free]");
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             SOLVER
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void compute_faces(mesh * inm, int ** edges, int sizee) {
    /*
     * This computes the half-edge data structure
     */
    const char *FUN = "compute_faces()" + supress_compute;
    inm->nofaces = inm->novertices + inm->notriangles;
    face * faces = (face *) malloc(sizeof(face) * inm->nofaces);
    int i, j, p, t;
    int tc[16] = {0};

    p = 0;
    t = -1;
    inm->noborders = 0;

    fori (sizee) {
        if (edges[i][2] != t) {
            inm->noborders++;
            t = edges[i][2];
        }
        tc[edges[i][2]]++;

        forj (inm->notriangles) {
            if (edges[i][0] == indexof(inm->vertices, inm->triangles[j].a) &&
                edges[i][1] == indexof(inm->vertices, inm->triangles[j].b)) {
                faces[p].vi = inm->vertices + edges[i][0];
                faces[p].vf = inm->vertices + edges[i][1];
                faces[p].tl = inm->triangles + j;
                faces[p].tr = 0;
                faces[p].vl = inm->triangles[j].c;
                faces[p].vr = 0;
                faces[p].border = edges[i][2];
                faces[p].vi->border = edges[i][2];
                faces[p].vf->border = edges[i][2];
                faces[p].vi->posface = faces + p;
                faces[p].vf->preface = faces + p;
                inm->triangles[j].border = edges[i][2];
                inm->triangles[j].curved = AB;
                inm->triangles[j].curvedface = faces + p;
                p++;
                break;
            }
            if (edges[i][0] == indexof(inm->vertices, inm->triangles[j].b) &&
                edges[i][1] == indexof(inm->vertices, inm->triangles[j].c)) {
                faces[p].vi = inm->vertices + edges[i][0];
                faces[p].vf = inm->vertices + edges[i][1];
                faces[p].tl = inm->triangles + j;
                faces[p].tr = 0;
                faces[p].vl = inm->triangles[j].a;
                faces[p].vr = 0;
                faces[p].border = edges[i][2];
                faces[p].vi->border = edges[i][2];
                faces[p].vf->border = edges[i][2];
                faces[p].vi->posface = faces + p;
                faces[p].vf->preface = faces + p;
                inm->triangles[j].border = edges[i][2];
                inm->triangles[j].curved = BC;
                inm->triangles[j].curvedface = faces + p;
                p++;
                break;
            }
            if (edges[i][0] == indexof(inm->vertices, inm->triangles[j].c) &&
                edges[i][1] == indexof(inm->vertices, inm->triangles[j].a)) {
                faces[p].vi = inm->vertices + edges[i][0];
                faces[p].vf = inm->vertices + edges[i][1];
                faces[p].tl = inm->triangles + j;
                faces[p].tr = 0;
                faces[p].vl = inm->triangles[j].b;
                faces[p].vr = 0;
                faces[p].border = edges[i][2];
                faces[p].vi->border = edges[i][2];
                faces[p].vf->border = edges[i][2];
                faces[p].vi->posface = faces + p;
                faces[p].vf->preface = faces + p;
                inm->triangles[j].border = edges[i][2];
                inm->triangles[j].curved = CA;
                inm->triangles[j].curvedface = faces + p;
                p++;
                break;
            }
        }
    }

    fori (inm->notriangles) {
        forjn (i +1, inm->notriangles) {
            if (inm->triangles[i].a == inm->triangles[j].b &&
                inm->triangles[i].b == inm->triangles[j].a) {
                faces[p].vi = inm->triangles[i].a;
                faces[p].vf = inm->triangles[i].b;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].c;
                faces[p].vr = inm->triangles[j].c;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].a == inm->triangles[j].c &&
                inm->triangles[i].b == inm->triangles[j].b) {
                faces[p].vi = inm->triangles[i].a;
                faces[p].vf = inm->triangles[i].b;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].c;
                faces[p].vr = inm->triangles[j].a;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].a == inm->triangles[j].a &&
                inm->triangles[i].b == inm->triangles[j].c) {
                faces[p].vi = inm->triangles[i].a;
                faces[p].vf = inm->triangles[i].b;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].c;
                faces[p].vr = inm->triangles[j].a;
                faces[p].border = NONE;
                p++;
                break;
            }
        }
        forjn (i +1, inm->notriangles) {
            if (inm->triangles[i].b == inm->triangles[j].b &&
                inm->triangles[i].c == inm->triangles[j].a) {
                faces[p].vi = inm->triangles[i].b;
                faces[p].vf = inm->triangles[i].c;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].a;
                faces[p].vr = inm->triangles[j].c;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].b == inm->triangles[j].c &&
                inm->triangles[i].c == inm->triangles[j].b) {
                faces[p].vi = inm->triangles[i].b;
                faces[p].vf = inm->triangles[i].c;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].a;
                faces[p].vr = inm->triangles[j].a;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].b == inm->triangles[j].a &&
                inm->triangles[i].c == inm->triangles[j].c) {
                faces[p].vi = inm->triangles[i].b;
                faces[p].vf = inm->triangles[i].c;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].a;
                faces[p].vr = inm->triangles[j].b;
                faces[p].border = NONE;
                p++;
                break;
            }
        }
        forjn (i +1, inm->notriangles) {
            if (inm->triangles[i].c == inm->triangles[j].b &&
                inm->triangles[i].a == inm->triangles[j].a) {
                faces[p].vi = inm->triangles[i].c;
                faces[p].vf = inm->triangles[i].a;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].b;
                faces[p].vr = inm->triangles[j].c;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].c == inm->triangles[j].c &&
                inm->triangles[i].a == inm->triangles[j].b) {
                faces[p].vi = inm->triangles[i].c;
                faces[p].vf = inm->triangles[i].a;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].b;
                faces[p].vr = inm->triangles[j].a;
                faces[p].border = NONE;
                p++;
                break;
            }
            if (inm->triangles[i].c == inm->triangles[j].a &&
                inm->triangles[i].a == inm->triangles[j].c) {
                faces[p].vi = inm->triangles[i].c;
                faces[p].vf = inm->triangles[i].a;
                faces[p].tl = inm->triangles + i;
                faces[p].tr = inm->triangles + j;
                faces[p].vl = inm->triangles[i].b;
                faces[p].vr = inm->triangles[j].b;
                faces[p].border = NONE;
                p++;
                break;
            }
        }
    }

    inm->faces = faces;

    inm->noedgeface = cint(inm->noborders);
    inm->edgeface   = (face **) malloc(sizeof(face *) * inm->noborders);

    inm->edgeface[0] = faces;
    fori (inm->noborders) {
        inm->noedgeface[i] = tc[i];
        if (i > 0)
            inm->edgeface[i] = faces + tc[i - 1];
    }

    if (!inm->quiet) INFOMF("proccessed %d faces (%d edges; %d different borders)", p, sizee, inm->noborders);
}

void compute_radius(mesh * inm) {
    /*
     * This computes the curvature radius of each edge face
     */
    const char *FUN = "compute_radius()" + supress_compute;
     int i, j, t = 0;

    /* for each border, compute the spline and more... */
    fori (inm->noborders) {
        double *tt = mdouble(inm->noedgeface[i] + 1),
               *xx = mdouble(inm->noedgeface[i] + 1), /* xx values must be strictly increasing */
               *yy = mdouble(inm->noedgeface[i] + 1),
               *curve_radius;

        tt[0] = 0;
        xx[0] = inm->edgeface[i][0].vi->x;
        yy[0] = inm->edgeface[i][0].vi->y;
        forjn (1, inm->noedgeface[i] + 1) {
            tt[j] = j;
            xx[j] = inm->edgeface[i][j - 1].vf->x;
            yy[j] = inm->edgeface[i][j - 1].vf->y;
        }

        curve_radius = compute_radius_spline(inm->noedgeface[i] + 1, tt, xx, yy) ;

        inm->edgeface[i][0].curve_radius                    = curve_radius[1];                    /* first radius */
        inm->edgeface[i][inm->noedgeface[i]-1].curve_radius = curve_radius[inm->noedgeface[i]-1]; /* last  radius */

        forjn (1, inm->noedgeface[i] - 1) {
            inm->edgeface[i][j].curve_radius = 0.5 * (curve_radius[j] + curve_radius[j+1]);
        }
        t += inm->noedgeface[i];

        if (compute_wall_condiction(inm->edgeface[i][0].border, inm->inflow_angle, 0.0, 0.0, 1.0, 0.0) == WALL) { /* WALL CONDICTION */
            forj (inm->noedgeface[i] - 1) {
                inm->edgeface[i][j].vf->curve_radius = curve_radius[j+1];
                inm->edgeface[i][j].vf->flow = WALL;
            }
        }

        /*printf("\n\n");
        forj (inm->noedgeface[i]) {
            printf("%f (vi: %f; vf: %f)\n", inm->edgeface[i][j].curve_radius, inm->edgeface[i][j].vi->curve_radius, inm->edgeface[i][j].vf->curve_radius);
        }
        printf("\n\n");*/

        /* Free! */
        free(tt);
        free(xx);
        free(yy);
        free(curve_radius);
    }
    if (!inm->quiet && inm->verbosity)
        INFOMF("%d splines computed as well as %d radii", inm->noborders, t);
}

void compute_middle_points(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_middle_points()" + supress_compute;*/
    int i;
    fori (inm->nofaces) {
        if (inm->faces[i].border == NONE) { /* NON-BORDER */
            inm->faces[i].middle_point.x = 0.5 * (inm->faces[i].vi->x + inm->faces[i].vf->x);
            inm->faces[i].middle_point.y = 0.5 * (inm->faces[i].vi->y + inm->faces[i].vf->y);
        } else { /* BORDER FACE */
            double tx = inm->faces[i].vf->x - inm->faces[i].vi->x;
            double ty = inm->faces[i].vf->y - inm->faces[i].vi->y;
            double l  = sqrt(tx * tx + ty * ty) ;
            tx = tx / l ;
            ty = ty / l ;

            double theta = asin (0.5 * l / inm->faces[i].curve_radius);
            inm->faces[i].middle_point.x = inm->faces[i].vf->x +
                    inm->faces[i].curve_radius * ((1 - cos(theta)) *  ty - sin(theta) * tx);
            inm->faces[i].middle_point.y = inm->faces[i].vf->y +
                    inm->faces[i].curve_radius * ((1 - cos(theta)) * -tx - sin(theta) * ty);
        }
    }
}

void compute_legendre(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_legendre()" + supress_compute;*/
    int n;

    inm->gauss_points [0] = mdouble(1);
    inm->gauss_weights[0] = mdouble(1);

    inm->gauss_points [0][0] = 0.0;
    inm->gauss_weights[0][0] = 2.0;

    fornn (2, 4) {
        int i , j , k , pos , q , r ;
        double p0[4] = {0}, p1[4] = {0}, a[4] = {0}, z[6] = {0}, vaux[4] = {0},
               an,  bn, min;

        zs(i, p0, 4); zs(i, p1, 4); zs(i, a, 4); zs(i, z, 6); zs(i, vaux, 4);

        inm->gauss_points [n - 1] = mdouble(n);
        inm->gauss_weights[n - 1] = mdouble(n);

        /* Legendre's polinomial */
        p0 [ 0 ] = 1.0 ;
        p1 [ 1 ] = 1.0 ;

        for (i = 1; i <= n - 1; i++) {
            an = -1.0 * i / (i + 1);
            bn = (2.0 * i + 1.0) / (i + 1);
            for (j = 0; j <= i; j++) {
                vaux[j + 1] = p1[j];
            }
            for (j = 0 ; j <= i + 1; j++) {
                a[j] = an * p0[j] + bn * vaux[j];
            }
            if (i < n - 1) {
                for (j = 0; j <= i + 1; j++) {
                    p0[j] = p1[j];
                    p1[j] = a[j];
                    a[j]  = 0.0;
                }
            }
        }

        /* Legendre's roots */
        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (n + 1) ;
        gsl_poly_complex_solve(a, n + 1, w, z);
        gsl_poly_complex_workspace_free(w) ;

        k = -1;

        /* Legendre's points */
        for (i = 0; i < n; i++) {
            ++k ; // k=? why not i?
            inm->gauss_points[n - 1][k] = z[2 * i];
        }

        for (i = 0; i < n - 1; i++) {
            pos = i ;
            min = inm->gauss_points[n - 1][i] ;
            for (j = i + 1; j < n; j++) {
                if (min > inm->gauss_points[n - 1][j]) {
                    pos = j;
                    min = inm->gauss_points[n - 1][j];
                }
            }
            if (pos > i) {
                min = inm->gauss_points[n - 1][pos];
                inm->gauss_points[n - 1][pos] = inm->gauss_points[n - 1][i];
                inm->gauss_points[n - 1][i]   = min;
            }
        }

        q = n / 2;
        r = n - 2 * q;

        if (r == 0) {
            for (i = 0; i < q; i++) {
                inm->gauss_points[n - 1][n - i - 1] = -inm->gauss_points[n - 1][i];
            }
        }
        else {
            for (i = 0; i < q; i++) {
                inm->gauss_points[n - 1][n - i - 1] = -inm->gauss_points[n - 1][i];
            }
            inm->gauss_points[n - 1][(n - 1) / 2] = 0.0;
        }

        /* Legendre's points (again?) */
        for (i = 0; i <= n; i++)
            a[i] = i * a[i];

        for (i = 1; i <= n; i++)
            a[i - 1] = a[i] ;

        a[n] = 0.0 ;

        for (i = 0; i < n; i++) {
            bn = a[n - 1];
            for (j = n - 2; j >= 0; j--) {
                bn = a[j] + bn * inm->gauss_points[n - 1][i];
            }
            /* weights */
            inm->gauss_weights[n - 1][i] = 2.0 / ((1.0 - inm->gauss_points[n - 1][i] * inm->gauss_points[n - 1][i]) * pow(bn, 2.0));
        }
    }
}

void compute_centroids(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_centroids()" + supress_compute;*/
    int i, j, n;

    fori(inm->notriangles) { /* for each triangle */
        if (inm->triangles[i].curved == NN) { /* INTERNAL TRIANGLE == no curved faces */
            inm->triangles[i].centroid.x =
                (inm->triangles[i].a->x + inm->triangles[i].b->x + inm->triangles[i].c->x) / 3.0;
            inm->triangles[i].centroid.y =
                (inm->triangles[i].a->y + inm->triangles[i].b->y + inm->triangles[i].c->y) / 3.0;
        }
        else {; /* BORDER TRIANGLE == one curved face */
            double momentum[3] = {0};
            forjn (AB, CA + 1) {
                double xa, ya, xb, yb;
                /* get the vertices' values of the face */
                switch (j) {
                    case AB:
                        xa = inm->triangles[i].a->x; ya = inm->triangles[i].a->y;
                        xb = inm->triangles[i].b->x; yb = inm->triangles[i].b->y;
                        break;
                    case BC:
                        xa = inm->triangles[i].b->x; ya = inm->triangles[i].b->y;
                        xb = inm->triangles[i].c->x; yb = inm->triangles[i].c->y;
                        break;
                    case CA:
                        xa = inm->triangles[i].c->x; ya = inm->triangles[i].c->y;
                        xb = inm->triangles[i].a->x; yb = inm->triangles[i].a->y;
                        break;
                }

                if (inm->triangles[i].curved != j) { /* NON-curved face integral computation */
                    forn (2) { /* two gauss points */
                        double xxgg = xa + 0.5 * (xb - xa) * (inm->gauss_points[1][n] + 1.0);
                        double yygg = ya + 0.5 * (yb - ya) * (inm->gauss_points[1][n] + 1.0);
                        double aux1 = inm->gauss_weights[1][n] * 0.5 * (yb - ya);
                        double aux2 = xxgg - inm->triangles[i].a->x;
                        double aux3 = yygg - inm->triangles[i].a->y;
                        momentum[0] += aux1 * aux2;
                        momentum[1] += aux1 * 0.5 * aux2 * aux2;
                        momentum[2] += aux1 * aux2 * aux3;
                    }
                } else { /* curved face integral computation */
                    double tx = xb - xa;
                    double ty = yb - ya;
                    double l  = sqrt(tx * tx + ty * ty) ;
                    tx = tx / l;
                    ty = ty / l;
                    double theta = asin(0.5 * l / inm->triangles[i].curvedface->curve_radius);

                    forn (3) { /* three gauss points */
                        double beta = theta * inm->gauss_points[2][n];
                        double xxgg = xb + inm->triangles[i].curvedface->curve_radius * ((cos(beta) - cos(theta)) *  ty + (sin(beta) - sin(theta)) * tx);
                        double yygg = yb + inm->triangles[i].curvedface->curve_radius * ((cos(beta) - cos(theta)) * -tx + (sin(beta) - sin(theta)) * ty);
                        double nxgg = cos(beta) * ty + sin(beta) * tx;
                        double wwgg = inm->gauss_weights[2][n] * fabs(inm->triangles[i].curvedface->curve_radius * theta); /* TODO abs with fabs */
                        double aux1 = wwgg * nxgg;
                        double aux2 = xxgg - inm->triangles[i].a->x;
                        double aux3 = yygg - inm->triangles[i].a->y;
                        momentum[0] += aux1 * aux2;
                        momentum[1] += aux1 * 0.5 * aux2 * aux2;
                        momentum[2] += aux1 * aux2 * aux3;
                    }
                }
            }
            inm->triangles[i].centroid.x = (momentum[1] / momentum[0]) + inm->triangles[i].a->x;
            inm->triangles[i].centroid.y = (momentum[2] / momentum[0]) + inm->triangles[i].a->y;
        }
    }
}


void compute_momentum(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_momentum()" + supress_compute;*/
    int i, j;

    /* allocte mem for momentum */
    fori (inm->novertices)
        inm->vertices[i].momentum = cdouble(inm->notc);

    fori (inm->nofaces) {
        if (inm->faces[i].border == NONE) { /* NON-border face */

            /* initial's face vertex */
            compute_momentum_vertex(inm->faces[i].vi->momentum, inm->gauss_points, inm->gauss_weights,
                    inm->order, inm->nommt, 0,
                    *((coord*) inm->faces[i].vi), inm->faces[i].tr->centroid, inm->faces[i].middle_point, 0);
            compute_momentum_vertex(inm->faces[i].vi->momentum, inm->gauss_points, inm->gauss_weights,
                    inm->order, inm->nommt, 0,
                    *((coord*) inm->faces[i].vi), inm->faces[i].middle_point, inm->faces[i].tl->centroid, 0);

            /* final's face vertex */
            compute_momentum_vertex(inm->faces[i].vf->momentum, inm->gauss_points, inm->gauss_weights,
                    inm->order, inm->nommt, 0,
                    *((coord*) inm->faces[i].vf), inm->faces[i].tl->centroid, inm->faces[i].middle_point, 0);
            compute_momentum_vertex(inm->faces[i].vf->momentum, inm->gauss_points, inm->gauss_weights,
                    inm->order, inm->nommt, 0,
                    *((coord*) inm->faces[i].vf), inm->faces[i].middle_point, inm->faces[i].tr->centroid, 0);
        } else {                            /* BORDER face     */
            /* initial's face vertex */
            compute_momentum_vertex(inm->faces[i].vi->momentum, inm->gauss_points, inm->gauss_weights,
                                    inm->order, 3, 1,
                                    *((coord*) inm->faces[i].vi), *((coord*) inm->faces[i].vi), inm->faces[i].middle_point, inm->faces[i].curve_radius);
            compute_momentum_vertex(inm->faces[i].vi->momentum, inm->gauss_points, inm->gauss_weights,
                                    inm->order, inm->nommt, 0,
                                    *((coord*) inm->faces[i].vi), inm->faces[i].middle_point, inm->faces[i].tl->centroid, 0);

            /* final's face vertex */
            compute_momentum_vertex(inm->faces[i].vf->momentum, inm->gauss_points, inm->gauss_weights,
                                    inm->order, inm->nommt, 0,
                                    *((coord*) inm->faces[i].vf), inm->faces[i].tl->centroid, inm->faces[i].middle_point, 0);
            compute_momentum_vertex(inm->faces[i].vf->momentum, inm->gauss_points, inm->gauss_weights,
                                    inm->order, 3, 1,
                                    *((coord*) inm->faces[i].vf), inm->faces[i].middle_point, *((coord*) inm->faces[i].vf), inm->faces[i].curve_radius);
        }
    }

    fori (inm->novertices) {
        inm->vertices[i].cva = inm->vertices[i].momentum[0];
        inm->vertices[i].momentum[0] = 1.0;
        forjn (1, inm->notc)
            inm->vertices[i].momentum[j] /= inm->vertices[i].cva;
    }
}

void compute_projection(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_projection()" + supress_compute;*/
    int i;

    fori (inm->novertices) {
        inm->vertices[i].projection.x = 0.0;
        inm->vertices[i].projection.y = 0.0;
    }

    fori (inm->nofaces) {
        if (inm->faces[i].border == NONE) { /* NON-border face */
            double x1 = 0.5 * (inm->faces[i].vi->x + inm->faces[i].vf->x),
                   y1 = 0.5 * (inm->faces[i].vi->y + inm->faces[i].vf->y),
                   x2 = inm->faces[i].tl->centroid.x,
                   y2 = inm->faces[i].tl->centroid.y,
                   x3 = inm->faces[i].tr->centroid.x,
                   y3 = inm->faces[i].tr->centroid.y;

            double t  = fabs(y1 - y3) + fabs(y2 - y1);
            inm->faces[i].vi->projection.x += t;
            inm->faces[i].vf->projection.x += t;

            t  = fabs(x3 - x1) + fabs(x1 - x2);
            inm->faces[i].vi->projection.y += t;
            inm->faces[i].vf->projection.y += t;
        } else {                            /* BORDER face     */
           double x1 = inm->faces[i].middle_point.x,
                   y1 = inm->faces[i].middle_point.y,
                   x2 = inm->faces[i].tl->centroid.x,
                   y2 = inm->faces[i].tl->centroid.y,
                   x3 = inm->faces[i].vi->x,
                   y3 = inm->faces[i].vi->y;

            inm->faces[i].vi->projection.x += fabs(y1 - y3);
            inm->faces[i].vi->projection.y += fabs(x3 - x1);

            double t = fabs(y2 - y1);
            inm->faces[i].vi->projection.x += t;
            inm->faces[i].vf->projection.x += t;

            t = fabs (x1 - x2);
            inm->faces[i].vi->projection.y += t;
            inm->faces[i].vf->projection.y += t;

            x3 = inm->faces[i].vf->x,
            y3 = inm->faces[i].vf->y;
            inm->faces[i].vf->projection.x += fabs(y1 - y3);
            inm->faces[i].vf->projection.y += fabs(x3 - x1);
        }
    }

    fori (inm->novertices) {
        inm->vertices[i].projection.x *= 0.5;
        inm->vertices[i].projection.y *= 0.5;
    }
}

void compute_stencil(mesh * inm) {
    /*
     * This computes TODO
     */
    const char *FUN = "compute_stencil()" + supress_compute;
    int i;

    if (!inm->quiet && inm->verbosity) WARNMF("[NIY-build]");

    /*#pragma message "compute_stencil(): NOT FULLY IMPLEMENTED!"*/

    /* allocate space for the stencils. TODO is that need so much mem? at least find an upper bound */
    vertex ** t = (vertex **) malloc(sizeof(vertex *) * inm->novertices * inm->novertices);
    fori (inm->novertices)
        inm->vertices[i].cv_stencil = t + inm->novertices * i;

    fori (inm->nofaces) {
        inm->faces[i].vi->cv_stencil[inm->faces[i].vi->cv_no] = inm->faces[i].vf;
        inm->faces[i].vf->cv_stencil[inm->faces[i].vf->cv_no] = inm->faces[i].vi;
        inm->faces[i].vi->cv_no++;
        inm->faces[i].vf->cv_no++;
    }

    int *posold = cint(inm->novertices),
        *posnew = cint(inm->novertices);
    fori (inm->novertices) {
        posold[i] = 1;
        posnew[i] = inm->vertices[i].cv_no;
        /*inm->vertices[i].cv_no_1l = inm->vertices[i].cv_no;*/
    }

    int cv_no_min = (int) ceil((inm->notc - 1.0) * 1.5); /*refence: Ollivier-Gooch et al 2009 */
    fori (inm->novertices) {
        if (inm->vertices[i].cv_no < cv_no_min) {
            AERRMF("Dude, do you wanna compute fluid dynamics without a mesh? Go build a fine tune one! [NIY-too-few-vertices]");
            /*int *m = cint(inm->novertices);
            m[0] = 1;

            forj (inm->vertices[i].cv_no_1l)
                m[inm->vertices[i].cv_stencil[j]] = 1;
            free(m);*/
        }
    }

    free(posold);
    free(posnew);
}

void compute_gausses(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_gausses()" + supress_compute;*/
    int i, j;

    fori (inm->nofaces) {
        coord vi, vf;
        double l;

        if (inm->faces[i].border == NONE) { /* NON-border face */
            /* Allocate mem */
            inm->faces[i].gauss_points  = (coord *)  malloc(sizeof(coord)  * 2 * inm->lpd);
            inm->faces[i].gauss_weights = (double *) malloc(sizeof(double) * 2 * inm->lpd);
            inm->faces[i].normals       = (coord *)  malloc(sizeof(coord)  * 2);
            inm->faces[i].flow          = (flow *)   malloc(sizeof(flow)   * 1);

            /* internal faces do not have special flow condiction */
            inm->faces[i].flow[0] = INTERNAL;

            /* first group */
            vi.x = inm->faces[i].tr->centroid.x;
            vi.y = inm->faces[i].tr->centroid.y;
            vf.x = 0.5 * (inm->faces[i].vi->x + inm->faces[i].vf->x);
            vf.y = 0.5 * (inm->faces[i].vi->y + inm->faces[i].vf->y);
            l = sqrt(pow2(vf.x - vi.x) + pow2(vf.y - vi.y));

            forj (inm->lpd) {
                inm->faces[i].gauss_points[j].x = vi.x + 0.5 * (vf.x - vi.x) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_points[j].y = vi.y + 0.5 * (vf.y - vi.y) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_weights[j]  = inm->gauss_weights[inm->lpd - 1][j] * 0.5 * l;
            }
            inm->faces[i].normals[0].x =  (vf.y - vi.y) / l;
            inm->faces[i].normals[0].y = -(vf.x - vi.x) / l;

            /* second group */
            vi.x = vf.x;
            vi.y = vf.y;
            vf.x = inm->faces[i].tl->centroid.x;
            vf.y = inm->faces[i].tl->centroid.y;
            l = sqrt(pow2(vf.x - vi.x) + pow2(vf.y - vi.y));

            forj (inm->lpd) {
                inm->faces[i].gauss_points[inm->lpd + j].x = vi.x + 0.5 * (vf.x - vi.x) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_points[inm->lpd + j].y = vi.y + 0.5 * (vf.y - vi.y) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_weights[inm->lpd + j]  = inm->gauss_weights[inm->lpd - 1][j] * 0.5 * l;
            }
            inm->faces[i].normals[1].x =  (vf.y - vi.y) / l;
            inm->faces[i].normals[1].y = -(vf.x - vi.x) / l;
        } else {                            /* BORDER face     */
            /* Allocate mem */
            inm->faces[i].gauss_points  = (coord *)  malloc(sizeof(coord)  * 3 * inm->lpd);
            inm->faces[i].gauss_weights = (double *) malloc(sizeof(double) * 3 * inm->lpd);
            inm->faces[i].normals       = (coord *)  malloc(sizeof(coord)  * 3 * inm->lpd);
            inm->faces[i].flow          = (flow *)   malloc(sizeof(flow)   * 3 * inm->lpd);

            double tx, ty, radius, theta, beta;

            /* first */
            vi.x   = inm->faces[i].vi->x;
            vi.y   = inm->faces[i].vi->y;
            vf.x   = inm->faces[i].middle_point.x;
            vf.y   = inm->faces[i].middle_point.y;
            radius = inm->faces[i].curve_radius;

            tx = vf.x - vi.x ;
            ty = vf.y - vi.y ;
            l  = sqrt(pow2(tx) + pow2(ty));
            tx = tx / l;
            ty = ty / l;
            theta = asin(0.5 * l / radius);

            forj (inm->lpd) {
                beta = theta * inm->gauss_points[inm->lpd - 1][j];
                inm->faces[i].gauss_points[j].x = vf.x + radius * ((cos(beta) - cos(theta)) *  ty + (sin(beta) - sin(theta)) * tx);
                inm->faces[i].gauss_points[j].y = vf.y + radius * ((cos(beta) - cos(theta)) * -tx + (sin(beta) - sin(theta)) * ty);
                inm->faces[i].gauss_weights[j]  = inm->gauss_weights[inm->lpd - 1][j] * fabs(radius * theta);
                inm->faces[i].normals[j].x      = cos(beta) *  ty + sin(beta) * tx;
                inm->faces[i].normals[j].y      = cos(beta) * -tx + sin(beta) * ty;
            }

            /* second */
            vi.x = inm->faces[i].middle_point.x;
            vi.y = inm->faces[i].middle_point.y;
            vf.x = inm->faces[i].tl->centroid.x;
            vf.y = inm->faces[i].tl->centroid.y;
            l    = l2dist(vi.x - vf.x, vi.y - vf.y);
            forj (inm->lpd) {
                int jo = j + 1 * inm->lpd;
                inm->faces[i].gauss_points[jo].x = vi.x + 0.5 * (vf.x - vi.x) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_points[jo].y = vi.y + 0.5 * (vf.y - vi.y) * (inm->gauss_points[inm->lpd - 1][j] + 1.0);
                inm->faces[i].gauss_weights[jo]  = inm->gauss_weights[inm->lpd - 1][j] * 0.5 * l;
                inm->faces[i].normals[jo].x      =  (vf.y - vi.y) / l;
                inm->faces[i].normals[jo].y      = -(vf.x - vi.x) / l;
            }

            /* third */
            vi.x = inm->faces[i].middle_point.x;
            vi.y = inm->faces[i].middle_point.y;
            vf.x = inm->faces[i].vf->x;
            vf.y = inm->faces[i].vf->y;

            tx = vf.x - vi.x ;
            ty = vf.y - vi.y ;
            l  = sqrt(pow2(tx) + pow2(ty));
            tx = tx / l;
            ty = ty / l;
            theta = asin(0.5 * l / radius);

            forj (inm->lpd) {
                int jo = j + 2 * inm->lpd;
                beta = theta * inm->gauss_points[inm->lpd - 1][j];
                inm->faces[i].gauss_points[jo].x = vf.x + radius * ((cos(beta) - cos(theta)) *  ty + (sin(beta) - sin(theta)) * tx);
                inm->faces[i].gauss_points[jo].y = vf.y + radius * ((cos(beta) - cos(theta)) * -tx + (sin(beta) - sin(theta)) * ty);
                inm->faces[i].gauss_weights[jo]  = inm->gauss_weights[inm->lpd - 1][j] * fabs(radius * theta);
                inm->faces[i].normals[jo].x      = cos(beta) *  ty + sin(beta) * tx;
                inm->faces[i].normals[jo].y      = cos(beta) * -tx + sin(beta) * ty;
            }
            forj (inm->lpd) {
                int j2l = 2 * inm->lpd + j;
                inm->faces[i].flow[j]   = compute_wall_condiction(inm->faces[i].border, inm->inflow_angle, inm->faces[i].gauss_points[j].x, inm->faces[i].gauss_points[j].y, inm->faces[i].normals[j].x, inm->faces[i].normals[j].y);
                inm->faces[i].flow[j2l] = compute_wall_condiction(inm->faces[i].border, inm->inflow_angle, inm->faces[i].gauss_points[j2l].x, inm->faces[i].gauss_points[j2l].y, inm->faces[i].normals[j2l].x, inm->faces[i].normals[j2l].y);
            }
        }
    }

}

void compute_coefmat_pseudoinverse(mesh * inm) {
    /*
     * This computes TODO
     */
    const char *FUN = "compute_coefmat_pseudoinverse()" + supress_compute;
    int i, j, n;
    double **A[inm->novertices];

    inm->max_cvno = 0;

    fori (inm->novertices) {
        /* allocate mem */
        double *t = cdouble((inm->vertices[i].cv_no + 1) * inm->notc);
        A[i]      = mdoublep(inm->vertices[i].cv_no + 1);

        forj (inm->vertices[i].cv_no + 1) //{
            A[i][j] = t + j * inm->notc; //printf("alloccing A[%d][%d]\n", i, j);}

        /* A[0][*] = momentum[*] */
        forj (inm->notc) {
            A[i][0][j] = inm->vertices[i].momentum[j];
        }
    }

    int ii, jj;
    double L;
    for (ii = 0; ii < inm->novertices; ii++) {
        inm->max_cvno = inm->vertices[ii].cv_no > inm->max_cvno ? inm->vertices[ii].cv_no : inm->max_cvno;
        for (jj = 0; jj < inm->vertices[ii].cv_no; jj++) {
            vertex *e = inm->vertices[ii].cv_stencil[jj];

            L = l2dist(inm->vertices[ii].x - e->x, inm->vertices[ii].y - e->y);
            L = 1.0 / pow(L, inm->Aexp);
            A[ii][jj + 1][0] = L;

            int c = 0;
            double xy[inm->order][inm->order];
            fori (inm->order) {
                forj (i + 1) {
                    xy[i - j][j] = e->momentum[c++];
                }
            }

            forin (1, inm->order) {
                forj (i + 1) {
                    int n = i - j;
                    int m = j;
                    double sum = 0.0;

                    int r, s;
                    for (r = 0; r < m + 1; r++)
                        for (s = 0; s < n + 1; s++) {
                            double bin1 = gsl_sf_fact(m) / (gsl_sf_fact(m - r) * gsl_sf_fact(r));
                            double bin2 = gsl_sf_fact(n) / (gsl_sf_fact(n - s) * gsl_sf_fact(s));
                            double p1   = gsl_pow_int(e->x - inm->vertices[ii].x, s);
                            double p2   = gsl_pow_int(e->y - inm->vertices[ii].y, r);
                            sum += bin1 * bin2 * p1 * p2 * xy[n - s][m - r];
                        }
                    A[ii][jj+1][j+1] = L * sum;
                }
            }
        }
    }

    if (!inm->quiet && inm->verbosity)
        INFOMF("Yo dawg, loop the loop the loop the... Let's pseudoinverse this");

    double *m[inm->max_cvno + 1];
    forn (inm->novertices) {
        /* alocate mem */
        /* m vector */
        inm->vertices[n].m = cdouble(inm->vertices[n].cv_no + 1);

        /* m temp matrix */
        double *t = cdouble(inm->vertices[n].cv_no * (inm->notc - 1));
        forj (inm->vertices[n].cv_no)
            m[j] = t + j * (inm->notc - 1);

        /* gauss' m vector and m matrix of unconstrained least square problem */
        forin (1, inm->vertices[n].cv_no + 1) {
            inm->vertices[n].m[i] = A[n][i][0];

            forjn (1, inm->notc) {
                m[i-1][j-1] = A[n][i][j] - inm->vertices[n].m[i] * A[n][0][j];
            }
        }

        /* Freeing up A */
        free(A[n][0]);
        free(A[n]);

        int l = inm->vertices[n].cv_no,
            c = inm->notc - 1;

        /* compute the pseudo-inverse */
        /* intput: l, c, m; output: inm->vertices[n].pinv */

        /* allocate mem for the pseudoinverse */
        t = cdouble(l * c);
        inm->vertices[n].pinv = mdoublep(l);
        forj (l)
            inm->vertices[n].pinv[j] = t + j * l;

        gsl_matrix_view  ma = gsl_matrix_view_array (m[0], l, c);
        gsl_matrix      *V  = gsl_matrix_alloc(c, c);
        gsl_vector      *S  = gsl_vector_alloc(c);
        gsl_vector      *w  = gsl_vector_alloc(c);
        gsl_linalg_SV_decomp (&ma.matrix, V, S, w);

        /* Nulling singular value if less than FEPS */
        const double FEPS = 1e-15;
        forj (c) {
            if (gsl_vector_get(S, j) < FEPS) {
                gsl_vector_set(S, j, 0.0);
            }
        }

        fori (l) {
            double *base = cdouble(l);
            base[i] = 1.0 ;
            gsl_vector_view  vecb = gsl_vector_view_array(base , l);
            gsl_vector      *xx   = gsl_vector_alloc(c);
            gsl_linalg_SV_solve(&ma.matrix, V, S, &vecb.vector, xx);

            forj (c) {
                inm->vertices[n].pinv[j][i] = gsl_vector_get(xx , j);
            }

            free(base);
            gsl_vector_free (xx) ;

        }
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_vector_free(w);

        /* Freeing up m temp */
        free(m[0]);
    }
}

void compute_prelimiters(mesh * inm) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_prelimiters()" + supress_compute;*/
    int i, j;

    /* Set auxiliary var for ghost pressure calc */
    fori (inm->novertices) {
        if (inm->vertices[i].flow == WALL) { /* A wall condiction vertex  */
            //printf("[%d] %e %e\n", i+1, inm->vertices[i].x, inm->vertices[i].y);

            coord vc_centroid = { inm->vertices[i].x + inm->vertices[i].momentum[1],
                                  inm->vertices[i].y + inm->vertices[i].momentum[2]};

            double d = l2dist(inm->vertices[i].momentum[1], inm->vertices[i].momentum[2]);
            forj (inm->lpd) {
                double deltax = inm->vertices[i].preface->gauss_points[2 * inm->lpd + j].x - vc_centroid.x;
                double deltay = inm->vertices[i].preface->gauss_points[2 * inm->lpd + j].y - vc_centroid.y;
                double aux    = l2dist(deltax, deltay);
                if (aux < d) d = aux;

                deltax = inm->vertices[i].posface->gauss_points[j].x - vc_centroid.x;
                deltay = inm->vertices[i].posface->gauss_points[j].y - vc_centroid.y;
                aux    = l2dist(deltax, deltay);
                if (aux < d) d = aux;
            }
            inm->vertices[i].a_ghostp = 2.0 * d / inm->vertices[i].curve_radius;
        }
    }

    /* Set auxiliary var: delta */
    double aux       = pow3(inm->k_delta);
    double aux_densi = pow2(inm->u_inflow[0]) * aux;
    double aux_veloc = inm->gamma * inm->r * inm->t_inflow * aux;
    double aux_press = pow2(inm->gamma * inm->u_inflow[3]) * aux;

    forj (inm->novertices) {
        aux = l2dist(inm->vertices[j].projection.x, inm->vertices[j].projection.y);
        aux = pow3(aux);
        inm->vertices[j].a_delta[0] = aux_densi * aux;
        inm->vertices[j].a_delta[1] = aux_veloc * aux;
        inm->vertices[j].a_delta[2] = inm->vertices[j].a_delta[1];
        inm->vertices[j].a_delta[3] = aux_press * aux;
    }

    /* find adjacent faces of each vertex */
    fori (inm->novertices) {
        if (inm->vertices[i].border == NONE) inm->vertices[i].aino = inm->vertices[i].cv_no;
        else                                 inm->vertices[i].aino = inm->vertices[i].cv_no - 2;
        inm->vertices[i].aiface = (face **) malloc(sizeof(face *) * inm->vertices[i].aino);
    }
    int *count = cint(inm->novertices);
    fori (inm->nofaces) {
        if (inm->faces[i].border != NONE)
            continue;

        int ii = indexof(inm->vertices, inm->faces[i].vi),
            fi = indexof(inm->vertices, inm->faces[i].vf);

        inm->faces[i].vi->aiface[count[ii]++] = inm->faces + i;
        inm->faces[i].vf->aiface[count[fi]++] = inm->faces + i;
    }
    free(count);
}

void compute_rungekutta5(mesh * inm) {
    /*
     * This computes TODO
     */
    const char *FUN = "compute_rungekutta5()" + supress_compute;
    int i, j, n;

    /* allocate mem for the uconserv */
    inm->uconserv       = mdoublep(inm->novertices);
    inm->uconserv[0]    = cdouble(inm->novertices * NOU);
    forin (1, inm->novertices)
        inm->uconserv[i] = inm->uconserv[0] + i * NOU;

    /* allocate mem for the u */
    inm->u              = mdoublep(inm->novertices);
    inm->u[0]           = cdouble(inm->novertices * NOU);
    forin (1, inm->novertices)
        inm->u[i] = inm->u[0] + i * NOU;

    /* allocate mem for the phi_chapel */
    inm->phi_chapel      = mdoublep(inm->novertices);
    inm->phi_chapel[0]   = cdouble(inm->novertices * NOU);
    forin (1, inm->novertices)
        inm->phi_chapel[i] = inm->phi_chapel[0] + i * NOU;

    /* alloc mem for coefs */
    inm->coef = (double ***) malloc(sizeof(double **) * inm->novertices);
    fori (inm->novertices) {
        inm->coef[i]    = mdoublep(NOU);
        inm->coef[i][0] = cdouble(NOU * inm->notc);
        forjn (1, NOU)
            inm->coef[i][j] = inm->coef[i][0] + j * inm->notc;
    }

    /* allocate mem for the r */
    inm->res       = mdoublep(inm->novertices);
    inm->res[0]    = cdouble(inm->novertices * NOU);
    forin (1, inm->novertices)
        inm->res[i] = inm->res[0] + i * NOU;

    /* helper pointers */
    double **uconserv   = inm->uconserv;
    double **u          = inm->u;
    double **phi_chapel = inm->phi_chapel;
    double ***coef      = inm->coef;
    double **r          = inm->res;

    /* temp mem */
    double dt[inm->novertices];
    double uconserv_0[inm->novertices][NOU];

    /* compute initial solution */
    fori (inm->novertices) {
        uconserv[i][0] = inm->u_inflow[0];
        uconserv[i][1] = inm->u_inflow[0] * inm->u_inflow[1];
        uconserv[i][2] = inm->u_inflow[0] * inm->u_inflow[2];
        uconserv[i][3] = inm->u_inflow[0] *
            (inm->r / (inm->gamma - 1.0) * inm->t_inflow + 0.5 * inm->q_inflow * inm->q_inflow);
    }

    /* init simulation - Runge-Kutta 5 stage : refer to Blazek 6.1.1 */
    int     iteration = 0;
    double  residue   = 1.0;
    double  rk_a[]    = { 0.0695, 0.1602, 0.2898, 0.5060};

    int maxs = (int) log10(inm->max_iter) + 1;

    //if (!inm->quiet && inm->showinner) printf("\n");
    running_solver = 1;
    while (residue > inm->nrt && iteration < inm->max_iter + 1 && !interrupt_solver) {
        //memcpy(uconserv_0[0], u[0], sizeof(double) * inm->novertices * NOU);
        fori (inm->novertices)
            forj (NOU)
                uconserv_0[i][j] = uconserv[i][j];

        compute_rk_convert(inm, uconserv, u); 
        compute_rk_reconstruction(inm, u, coef); 
        compute_rk_limiters(inm, u, coef, phi_chapel); 
        compute_rk_residue(inm, u, coef, phi_chapel, r, iteration); 

        /* set time step */
        fori (inm->novertices) {
            double c = sqrt(inm->gamma * u[i][3] / u[i][0]);
            double lambdaX = (fabs(u[i][1]) + c) * inm->vertices[i].projection.x;
            double lambdaY = (fabs(u[i][2]) + c) * inm->vertices[i].projection.y;
            dt[i] = inm->cfl * inm->vertices[i].cva / (lambdaX + lambdaY);
        }

        /* compute the residue */
        residue = 0.0;
        fori (inm->novertices) {
            residue += pow2(r[i][0]);
        }
        residue = sqrt(residue / inm->novertices);

        /* 5 stages of RK5 */
        forn (sizeof(rk_a) / sizeof(double)) {
            fori (inm->novertices)
                forj (NOU)
                    uconserv[i][j] = uconserv_0[i][j] - rk_a[n] * dt[i] * r[i][j];
            compute_rk_convert(inm, uconserv, u); //pr(u, inm->novertices, NOU); exit(0);
            compute_rk_reconstruction(inm, u, coef);
            compute_rk_limiters(inm, u, coef, phi_chapel);
            compute_rk_residue(inm, u, coef, phi_chapel, r, iteration);
        }
        fori (inm->novertices)
            forj (NOU)
                uconserv[i][j] = uconserv_0[i][j] - dt[i] * r[i][j];

        if (!inm->quiet && inm->showdetails && iteration % inm->showdetails == 0) {
            INFOMF("[%*d]: |r| = %.*e, Cd = %.*e, Cl = %.*e", maxs, iteration, ITRESSIZE, residue, ITRESSIZE, 0.0, ITRESSIZE, 0.0);
        }

        if (inm->showgraphics && ((iteration % inm->showgraphics == 0) || iteration == 0)) {
            v_draw_coefs(inm, 0);
        }

        iteration++;
    }
    running_solver = 0;
    iteration--;

    /*print2d("r", r, inm->novertices, NOU, 1);
    print2d("uconserv_0", uconserv_0, inm->novertices, NOU, 0);
    print2d("u", u, inm->novertices, NOU, 0);
    print2d("uconserv", uconserv, inm->novertices, NOU, 0);
    print2d("phi_chapel", phi_chapel, inm->novertices, NOU, 0);
    print3d(coef, inm->novertices, NOU, inm->notc, 0);*/


    if (!inm->quiet || inm->fclassify) {
        if (residue <= inm->nrt && iteration >= inm->max_iter + 1) {
            INFOMF("Finished after reaching the maximum iterations and reaching the residue threshold");
        } else if (residue <= inm->nrt) {
            INFOMF("Finished after reaching the residue threshold");
        } else if (iteration >= inm->max_iter + 1) {
            INFOMF("Finished after reaching the maximum iterations");
        } else if (residue != residue) {
            WARNMF("Finished after creating a black hole");
        } else { /* user interrupt */
            WARNMF("Finished after a ghost sent me a SIGINT!");
        }

        INFOMF("Finish after %d iteration%s", iteration, iteration > 1 ? "s" : "");
        INFOMF("Residue: %.*e", FRESSIZE, residue);
        WARNMF("Quality: [NIY-RK] Cp and Cl");

    }

    /* last update of grapgics */
    if (inm->showgraphics) {
        v_draw_coefs(inm, 0);
    }
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             HELPERS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
static void compute_rk_convert(mesh * inm, double **uc, double **u) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_convert()" + supress_compute;*/
    int i;
    fori (inm->novertices) {
        u[i][0] = uc[i][0];
        u[i][1] = uc[i][1] / uc [i][0];
        u[i][2] = uc[i][2] / uc [i][0];
        u[i][3] = (inm->gamma - 1.0) * (uc[i][3] - 0.5 * (uc[i][1] * uc[i][1] + uc[i][2] * uc[i][2]) / uc[i][0]);
    }
}

static void compute_rk_reconstruction(mesh * inm, double **u, double ***coef) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_reconstruction()" + supress_compute;*/
    int i, j, n, m;

    double vb[inm->novertices][NOU][inm->max_cvno + 1];

    fori (inm->novertices) {
        forj (NOU) {
            vb[i][j][0] = u[i][j];
        }
    }
    /* TODO merge this two fors */
    fori (inm->novertices) {
        forj (inm->vertices[i].cv_no) {
            vertex *s = inm->vertices[i].cv_stencil[j]; // k
            double L = 1.0 / pow(l2dist(inm->vertices[i].x - s->x, inm->vertices[i].y - s->y), inm->Aexp);
            forn (NOU)
                vb[i][n][j + 1] = L * u[indexof(inm->vertices, s)][n];
        }
    }

    double sum,
           bb[inm->max_cvno],
           xx[inm->notc - 1];

    forn (inm->novertices) {
        form (NOU) {

            forin (1, inm->vertices[n].cv_no + 1) {
                bb [i - 1] = vb[n][m][i] - inm->vertices[n].m[i] * vb[n][m][0];
            }

            /* xx = multply pseuvoinv by vector bb */
            fori (inm->notc - 1) {
                xx[i] = 0.0;
                forj (inm->vertices[n].cv_no)
                    xx[i] +=  inm->vertices[n].pinv[i][j] * bb[j];
            }

            sum = 0.0;
            forjn (1, inm->notc) {
                coef[n][m][j] = xx[j - 1];
                sum += coef[n][m][j] * inm->vertices[n].momentum[j];
            }
            coef[n][m][0] = u[n][m] - sum;
        }
    }
}

static double compute_rk_polinomial(mesh * inm, double ***coef, int p, int k, double x, double y) {
    /*
     * This computes TODO   \PXY()/
     */
    /*const char *FUN = "compute_rk_polinomial()" + supress_compute;*/
    int i , j , q , n , m ;
    double sum = 0.0 ;
    q = 0 ;
    fori (inm->order) {
        forj (i + 1) {
            n = i - j;
            m = j;
            sum += coef[p][k][q++] * gsl_pow_int(x - inm->vertices[p].x, n) * gsl_pow_int(y - inm->vertices[p].y, m);
        }
    }
    return sum;
}

static void compute_rk_limiters(mesh * inm, double **u, double ***coef, double **phi_chapel) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_limiters()" + supress_compute;*/
    int i, j, n, m;

    const double FEPS = 1e-6; /* avoiding low pressure or negative pressure (wuuuut?) */

    double nmach[inm->novertices];

    fori (inm->novertices) {
        double q_i = l2dist(u[i][1], u[i][2]);
        nmach[i] = q_i / sqrt(inm->gamma * u[i][3] / u[i][0]);

        if (inm->vertices[i].flow == WALL) { /* wall condiction vertex */
            /* TODO check precedences */
            double Ptotal_i = u[i][3] * pow(1.0 + (inm->gamma - 1.0) / 2.0 * nmach[i] * nmach[i], inm->gamma / (inm->gamma-1.0)); /* total pressure */
            vertex *w = inm->vertices + i;    //p = WallindexofV[i] ;
            w->ghostp[3] = u[i][3] + w->a_ghostp * u[i][0] * q_i * q_i; /* ghost pressure */
            if (w->ghostp[3] < FEPS)
                w->ghostp[3] = FEPS;
            double nmach_g = sqrt((pow(Ptotal_i / w->ghostp[3], (inm->gamma-1.0) / inm->gamma) -1.0) * 2.0 / (inm->gamma - 1.0)); /* ghost mach number */
            w->ghostp[0] = u[i][0] * pow(w->ghostp[3] / u[i][3], 1.0 / inm->gamma); /* ghost density */
            double q_g = nmach_g * sqrt(inm->gamma * w->ghostp[3] / w->ghostp[0]); /* velocity ghost modulo */
            w->ghostp[1] = (q_g / q_i) * u[i][1]; /* x ghost velocity */
            w->ghostp[2] = (q_g / q_i) * u[i][2]; /* y ghost velocity */
        }
    }

    fori (inm->novertices) {
        vertex *vi = inm->vertices + i;

        /* max mach in stencil */
        double max_mach = nmach[i];
        forj (vi->cv_no) {
            double t = nmach[indexof(inm->vertices, vi->cv_stencil[j])];
            max_mach = t > max_mach ? t : max_mach;
        }
        /* set phi_chapel if this CV is a stagnation point */
        if (max_mach <= inm->M1) {
            forj (NOU) {
                phi_chapel[i][j] = 1.0;
            }
        }
        else { /* when CV is NOT a stagnation point */
            double aux, sigma_chapel_i = 0;


            if (max_mach < inm->M2 && max_mach > inm->M1) {
                aux = (max_mach - inm->M1) / (inm->M2 - inm->M1);
                sigma_chapel_i = aux * aux *(2.0 * aux - 3.0) + 1.0;
            }
            else if (max_mach >= inm->M2)
                sigma_chapel_i = 0.0;

            forn (NOU) {
                /* define deltaUmax_i, deltaUmin_i and deltaU_i2 for primitive var n  */
                double deltaUmax_i = 0.0;
                double deltaUmin_i = 0.0;
                double deltaU_i2;
                double phi_i_n, phi_til_i_n = 0, sigma_til_i_n = 0;
                double xx, yy, aux1, aux2;

                if (inm->vertices[i].flow == WALL) { /* wall condiction vertex */
                    aux = inm->vertices[i].ghostp[n] - u[i][n];
                    if      (aux > deltaUmax_i) deltaUmax_i = aux;
                    else if (aux < deltaUmin_i) deltaUmin_i = aux;
                }
                forj (vi->cv_no) { /* TODO: WARNING! should be for the FIRST LAYER of stencil */
                                   /*                remember that BY NOW the first layer is the same size : stencil() */
                    vertex *vsj = vi->cv_stencil[j];
                    aux = u[indexof(inm->vertices, vsj)][n] - u[i][n];
                    if      (aux > deltaUmax_i) deltaUmax_i = aux ;
                    else if (aux < deltaUmin_i) deltaUmin_i = aux ;

                    if (vsj->flow == WALL) { /* a wall condiction cv vertex */
                        aux = vsj->ghostp[n] - u[i][n] ;
                        if      (aux > deltaUmax_i) deltaUmax_i = aux;
                        else if (aux < deltaUmin_i) deltaUmin_i = aux;
                    }
                }

                deltaU_i2 = pow2(deltaUmax_i - deltaUmin_i);
                if (deltaU_i2 <= vi->a_delta[n]) /* set phi_chapel if deltaU_i2 is too small */
                    phi_chapel[i][n] = 1.0;
                else {                           /* or set the sigma_til_i_n if not too small */
                    if (deltaU_i2 < 2.0 * vi->a_delta[n] && deltaU_i2 > vi->a_delta[n]) {
                        aux = deltaU_i2 / vi->a_delta[n] - 1.0;
                        sigma_til_i_n = aux * aux * (2.0 * aux - 3.0) + 1.0;
                    }
                    else if (deltaU_i2 >= 2.0 * vi->a_delta[n])
                        sigma_til_i_n = 0.0 ;

                    /* and set phi_i_n */
                    phi_i_n = 1.0;
                    forj (vi->aino) {
                        face* fj = vi->aiface[j];
                        form (2 * inm->lpd) {
                            xx = fj->gauss_points[m].x;
                            yy = fj->gauss_points[m].y;
                            aux = compute_rk_polinomial(inm, coef, i, n, xx, yy) - u[i][n];
                            if      (aux > 0.0) aux1 = deltaUmax_i / aux;
                            else if (aux < 0.0) aux1 = deltaUmin_i / aux;
                            else                aux1 = 1.5; /* a magic number ? */

                            if (aux1 < 1.5 )  aux2 = aux1 * (-4.0 / 27.0 * aux1 * aux1 + 1.0 ) ;
                            else              aux2 = 1.0;

                            if (aux2 < phi_i_n) phi_i_n = aux2 ;
                        }
                    }
                    if (vi->border != NONE) { /* a border vertex */
                        int index;
                        face* f;
                        forj (2) {
                            if (j == 0) {
                                index = inm->lpd;
                                f = vi->preface;
                            }
                            else {
                                index = 0;
                                f = vi->posface;
                            }

                            formn (index, index + 2 * inm->lpd) {
                                xx = f->gauss_points[m].x;
                                yy = f->gauss_points[m].y;
                                aux = compute_rk_polinomial(inm, coef, i, n, xx, yy) - u[i][n];
                                if      (aux > 0.0) aux1 = deltaUmax_i / aux;
                                else if (aux < 0.0) aux1 = deltaUmin_i / aux;
                                else                aux1 = 1.5; /* a magic number ? */

                                if (aux1 < 1.5 )  aux2 = aux1 * (-4.0 / 27.0 * aux1 * aux1 + 1.0 ) ;
                                else              aux2 = 1.0;

                                if (aux2 < phi_i_n) phi_i_n = aux2 ;
                            }
                        }
                    }

                    phi_til_i_n      = sigma_til_i_n + (1.0 - sigma_til_i_n) * phi_i_n;
                    phi_chapel[i][n] = sigma_chapel_i + (1.0 - sigma_chapel_i) * phi_til_i_n;
                }
            }
        }
    }
}

static double compute_rk_polinomial_lim(mesh * inm, double **u , double ***coef, double **phi_chapel, int p, int k, double x, double y) {
    /*
     * This computes TODO (\LimPXY/)
     */
    /*const char *FUN = "compute_rk_polinomial_lim()" + supress_compute;*/
    int i, j, q, n, m;
    double sum = 0.0;
    q = 1;
    forin (1, inm->order) {
        forj (i + 1) {
            n = i - j;
            m = j;
            sum += coef[p][k][q] * (gsl_pow_int(x - inm->vertices[p].x, n) * gsl_pow_int(y - inm->vertices[p].y, m) - inm->vertices[p].momentum[q]);
            ++q ;
        }
    }
    sum = u[p][k] + phi_chapel[p][k] * sum;
    return sum;
}

static void compute_rk_r_primitive_lim(mesh * inm, double *vp, double **u , double ***coef, double **phi_chapel, int p, double x, double y) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_r_primitive_lim()" + supress_compute;*/
    int k ;
    forvar (k, 0, NOU)
        vp[k] = compute_rk_polinomial_lim(inm, u, coef, phi_chapel, p, k, x, y);
}

static void compute_rk_r_subinflow(mesh * inm, double * f, double press, double nx, double ny) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_r_subinflow()" + supress_compute;*/
    double Mach, Temp, q, u, v, densi;
    double ental, aux;
    Mach = sqrt((pow(inm->pt_inflow / press, (inm->gamma - 1.0) / inm->gamma) - 1.0) * 2.0 / (inm->gamma-1.0));
    Temp = inm->tt_inflow / (1.0 + (inm->gamma - 1.0) / 2.0 * pow2(Mach));
    q = Mach * sqrt(inm->gamma * inm->r * Temp);
    u = q * cos(inm->inflow_angle);
    v = q * sin(inm->inflow_angle);
    densi = press/(inm->r * Temp);
    ental = 0.5 * (pow2(u) + pow2(v)) + press * inm->gamma / (densi*(inm->gamma - 1.0));
    aux = u * nx + v * ny;
    f[0] = densi * aux;
    f[1] = f[0] * u + nx * press;
    f[2] = f[0] * v + ny * press;
    f[3] = f[0] * ental;
}

static void compute_rk_r_suboutflow(mesh * inm, double * f, double * vp, double nx, double ny) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_r_suboutflow()" + supress_compute;*/
    double densi, u, v, press;
    double ental, aux;
    densi = vp[0];
    u = vp[1];
    v = vp[2];
    press = inm->u_inflow[3];
    ental = 0.5 * (pow2(u) + pow2(v)) + press * inm->gamma / (densi * (inm->gamma - 1.0));
    aux = u * nx + v * ny;
    f[0] = densi * aux;
    f[1] = f[0] * u + nx * press;
    f[2] = f[0] * v + ny * press;
    f[3] = f[0] * ental;
}

static void compute_rk_r_roe(mesh * inm, double *F, double *vp1, double *vp2, double nx, double ny) {
    /*
     * This computes TODO
     */
    /*const char *FUN = "compute_rk_r_roe()" + supress_compute;*/
    double Gamma = inm->gamma;

    double densi[2], u[2], v[2], press[2];
    densi[0] = vp1[0];   densi[1] = vp2[0];
    u[0]     = vp1[1];   u[1]     = vp2[1];
    v[0]     = vp1[2];   v[1]     = vp2[2];
    press[0] = vp1[3];   press[1] = vp2[3];

    int i, k ;
    double ental[2], fluxo[4][2], aux;

    fori (2) {
        ental[i] = 0.5 * (pow2(u[i]) + pow2(v[i])) + press[i] * Gamma / (densi[i] * (Gamma - 1.0));
        aux = nx * u[i] + ny*v[i];
        fluxo[0][i] = densi[i] * aux;
        fluxo[1][i] = densi[i] * u[i] * aux + nx * press[i];
        fluxo[2][i] = densi[i] * v[i] * aux + ny * press[i];
        fluxo[3][i] = densi[i] * ental[i] * aux;
    }

    double delta_densi = densi[1] - densi[0],
           delta_u     = u[1] - u[0],
           delta_v     = v[1] - v[0],
           delta_press = press[1] - press[0],
           delta_Vel_n = nx * delta_u + ny * delta_v;

           aux       = sqrt(densi[0]) + sqrt(densi[1]);
    double densi_til = sqrt(densi[0] * densi[1]),
           u_til     = (sqrt(densi[0]) * u[0] + sqrt(densi[1]) * u[1]) / aux,
           v_til     = (sqrt(densi[0]) * v[0] + sqrt(densi[1]) * v[1]) / aux,
           ental_til = (sqrt(densi[0]) * ental[0] + sqrt(densi[1]) * ental[1]) / aux,
           sqr_til   = pow2(u_til) + pow2(v_til),
           c_til     = sqrt((Gamma - 1.0) * (ental_til - 0.5 * sqr_til)),
           Vel_n_til = nx * u_til + ny * v_til;

    double lambda[4], alpha[4], autoveK[4][4];

    lambda[0] = Vel_n_til - c_til;
    lambda[1] = Vel_n_til;
    lambda[2] = Vel_n_til;
    lambda[3] = Vel_n_til + c_til;

    autoveK[0][0] = 1.0;                       autoveK[0][1] = 1.0;         autoveK[0][2] = 0.0;               autoveK[0][3] = 1.0;
    autoveK[1][0] = u_til-nx*c_til;            autoveK[1][1] = u_til;       autoveK[1][2] = ny;                autoveK[1][3] = u_til+nx*c_til;
    autoveK[2][0] = v_til-ny*c_til;            autoveK[2][1] = v_til;       autoveK[2][2] = -1.0*nx;           autoveK[2][3] = v_til+ny*c_til;
    autoveK[3][0] = ental_til-c_til*Vel_n_til; autoveK[3][1] = 0.5*sqr_til; autoveK[3][2] = u_til*ny-v_til*nx; autoveK[3][3] = ental_til+c_til*Vel_n_til;

    alpha[0] = (delta_press-densi_til*c_til*delta_Vel_n) / (2.0*c_til*c_til);
    alpha[1] = delta_densi - delta_press/(c_til*c_til);
    alpha[2] = densi_til * (ny*delta_u-nx*delta_v);
    alpha[3] = (delta_press+densi_til*c_til*delta_Vel_n) / (2.0*c_til*c_til);

    double fluxo_roe[4] = {0};

    forvar (k, 0, 4) {
        fori (4) {
            fluxo_roe[k] += fabs(lambda[i])*alpha[i]*autoveK[k][i];
        }
        fluxo_roe[k]=0.5*(fluxo[k][0]+fluxo[k][1])-0.5*fluxo_roe[k];
    }
    F[0] = fluxo_roe[0];
    F[1] = fluxo_roe[1];
    F[2] = fluxo_roe[2];
    F[3] = fluxo_roe[3];
}

static void compute_rk_residue(mesh * inm, double **u, double ***coef, double **phi_chapel, double **r, int iter) {
    const char *FUN = "compute_rk_residue()" + supress_compute;
    int i, j, n;
    double xg, yg, nx, ny, press, Mach, aux;

    #define paredeinv(f, press, nx, ny) f[0] = 0.0; f[1] = (nx)*(press); f[2] = (ny)*(press); fluxoinv[3] = 0.0

    double fluxoinv[NOU] = {0},
           vp1[NOU] = {0},
           vp2[NOU] = {0};

    double **ifinv;
    ifinv       = mdoublep(inm->novertices);
    ifinv[0]    = cdouble(inm->novertices * NOU);
    forin (1, inm->novertices)
        ifinv[i] = ifinv[0] + i * NOU;

    fori (inm->nofaces) {
        vertex *vi = inm->faces[i].vi,  /* vc1 */
               *vf = inm->faces[i].vf;  /* vc2 */

        if (inm->faces[i].border != NONE) { /* a border face*/
            /* integral computation of vi */
            forj (inm->lpd) {
                xg = inm->faces[i].gauss_points[j].x;
                yg = inm->faces[i].gauss_points[j].y;
                nx = inm->faces[i].normals[j].x;
                ny = inm->faces[i].normals[j].y;

                if (inm->faces[i].flow[j] == WALL) { //fprintf(f, "%d wall\n", i+1);
                    press = compute_rk_polinomial_lim(inm, u, coef, phi_chapel,  indexof(inm->vertices, vi), 3, xg, yg);
                    paredeinv(fluxoinv, press, nx, ny);
                } else
                if (inm->faces[i].flow[j] == INFLOW) { //fprintf(f, "%d inflow\n", i+1);
                    if (inm->mach_inflow > 1.0) {
                        /* TODO NEVER THE CASE, by now */
                        AERRMF("Dude, you'r trying to fly over 1 MACH! You sould buy a better aircraft or it will crash [NIY-inflow-vi]");
                    } else {
                        press = compute_rk_polinomial_lim(inm, u, coef, phi_chapel, indexof(inm->vertices, vi), 3, xg, yg);
                        compute_rk_r_subinflow(inm, fluxoinv, press, nx, ny);
                    }
                } else
                if (inm->faces[i].flow[j] == OUTFLOW) { //fprintf(f, "%d outflow\n", i+1);
                    compute_rk_r_primitive_lim(inm, vp1, u, coef, phi_chapel, indexof(inm->vertices, vi), xg, yg);
                    Mach = l2dist(vp1[1], vp1[2]) / sqrt(inm->gamma * vp1[3] / vp1[0]);
                    if (Mach > 1.0) {
                        /* TODO NEVER THE CASE, by now */
                        AERRMF("Dude, you'r trying to fly over 1 MACH! You sould buy a better aircraft or it will crash [NIY-outflow-vi]");
                    } else {
                        compute_rk_r_suboutflow(inm, fluxoinv, vp1, nx, ny);
                    }
                }
                forn (NOU) {
                    ifinv[indexof(inm->vertices, vi)][n] += fluxoinv[n] * inm->faces[i].gauss_weights[j];
                }
            }

            /* integral computation of vi and vf */
            forjn (1 * inm->lpd, 2 * inm->lpd) {
                xg = inm->faces[i].gauss_points[j].x;
                yg = inm->faces[i].gauss_points[j].y;
                nx = inm->faces[i].normals[j].x;
                ny = inm->faces[i].normals[j].y;

                compute_rk_r_primitive_lim(inm, vp1, u, coef, phi_chapel, indexof(inm->vertices, vi), xg, yg);
                compute_rk_r_primitive_lim(inm, vp2, u, coef, phi_chapel, indexof(inm->vertices, vf), xg, yg);
                compute_rk_r_roe(inm, fluxoinv, vp1, vp2, nx, ny);

                forn (NOU) {
                    aux = fluxoinv[n] * inm->faces[i].gauss_weights[j];
                    ifinv[indexof(inm->vertices, vi)][n] += aux ;
                    ifinv[indexof(inm->vertices, vf)][n] += -1.0 * aux;
                }
            }

            /* integral computation of vf */
            forjn (2 * inm->lpd, 3 * inm->lpd) {
                xg = inm->faces[i].gauss_points[j].x;
                yg = inm->faces[i].gauss_points[j].y;
                nx = inm->faces[i].normals[j].x;
                ny = inm->faces[i].normals[j].y;

                if (inm->faces[i].flow[j] == WALL) {
                    press = compute_rk_polinomial_lim(inm, u, coef, phi_chapel,  indexof(inm->vertices, vf), 3, xg, yg);
                    paredeinv(fluxoinv, press, nx, ny);
                } else
                if (inm->faces[i].flow[j] == INFLOW) {
                    if (inm->mach_inflow > 1.0) {
                        /* TODO NEVER THE CASE, by now */
                        AERRMF("Dude, you'r trying to fly over 1 MACH! You sould buy a better aircraft or it will crash [NIY-inflow-vf]");
                    } else {
                        press = compute_rk_polinomial_lim(inm, u, coef, phi_chapel,  indexof(inm->vertices, vf), 3, xg, yg);
                        compute_rk_r_subinflow(inm, fluxoinv, press, nx, ny);
                    }
                } else
                if (inm->faces[i].flow[j] == OUTFLOW) {
                    compute_rk_r_primitive_lim(inm, vp2, u, coef, phi_chapel, indexof(inm->vertices, vf), xg, yg);
                    Mach = l2dist(vp1[1], vp1[2]) / sqrt(inm->gamma * vp1[3] / vp1[0]);
                    if (Mach > 1.0) {
                        /* TODO NEVER THE CASE, by now */
                        AERRMF("Dude, you'r trying to fly over 1 MACH! You sould buy a better aircraft or it will crash [NIY-outflow-vf]");
                    } else {
                        compute_rk_r_suboutflow(inm, fluxoinv, vp2, nx, ny);
                    }
                }
                forn (NOU)
                    ifinv[indexof(inm->vertices, vf)][n] += fluxoinv [n] * inm->faces[i].gauss_weights[j];
            }
        } else {
            nx = inm->faces[i].normals[0].x;
            ny = inm->faces[i].normals[0].y;
            forjn (0, inm->lpd) {
                xg = inm->faces[i].gauss_points[j].x;
                yg = inm->faces[i].gauss_points[j].y;

                compute_rk_r_primitive_lim(inm, vp1, u, coef, phi_chapel, indexof(inm->vertices, vi), xg, yg);
                compute_rk_r_primitive_lim(inm, vp2, u, coef, phi_chapel, indexof(inm->vertices, vf), xg, yg);
                compute_rk_r_roe(inm, fluxoinv, vp1, vp2, nx, ny);

                forn (NOU) {
                    aux = fluxoinv[n] * inm->faces[i].gauss_weights[j];
                    ifinv[indexof(inm->vertices, vi)][n] += aux ;
                    ifinv[indexof(inm->vertices, vf)][n] += -1.0 * aux;
                }
            }

            nx = inm->faces[i].normals[1].x;
            ny = inm->faces[i].normals[1].y;
            forjn (inm->lpd, 2 * inm->lpd) {
                xg = inm->faces[i].gauss_points[j].x;
                yg = inm->faces[i].gauss_points[j].y;

                compute_rk_r_primitive_lim(inm, vp1, u, coef, phi_chapel, indexof(inm->vertices, vi), xg, yg);
                compute_rk_r_primitive_lim(inm, vp2, u, coef, phi_chapel, indexof(inm->vertices, vf), xg, yg);
                compute_rk_r_roe(inm, fluxoinv, vp1, vp2, nx, ny);

                forn (NOU) {
                    aux = fluxoinv[n] * inm->faces[i].gauss_weights[j];
                    ifinv[indexof(inm->vertices, vi)][n] += aux ;
                    ifinv[indexof(inm->vertices, vf)][n] += -1.0 * aux;
                }
            }
        }
    }

    fori (inm->novertices) {
        forj (NOU) {
            r[i][j] = ifinv[i][j] / inm->vertices[i].cva;
        }
    }

    free(ifinv[0]);
    free(ifinv);
}

void compute_momentum_vertex(double * m, double ** gp, double ** gw, int o, int nommt, int t, coord vr, coord vi, coord vf, double r) {
    /* m         : where to save the computed result
       gp        : gauss points
       gw        : gauss weights
       o         : which order
       nommt     : number of gauss points to use
       t         : 0 = line, 1 = curve
       vr, vi, vf: reference, initial and final coordinates
       r         : radius of curve
     */
     int i, j, n, c;

     double xxgg, yygg, nxgg, wwgg,
        tx, ty, l, theta, beta,
        aux1, aux2, aux3;

    if (t == 0 || t == 2) { /* a line */
        fori (nommt) {
            xxgg = vi.x + 0.5 * (vf.x - vi.x) * (gp[nommt - 1][i] + 1.0);
            yygg = vi.y + 0.5 * (vf.y - vi.y) * (gp[nommt - 1][i] + 1.0);
            aux1 = gw[nommt - 1][i] * 0.5 * (vf.y - vi.y);
            aux2 = xxgg - vr.x;
            aux3 = yygg - vr.y;
            c = 0 ;
            forj (o)
                forn (j +1)
                    m[c++] += aux1 * gsl_pow_int(aux2 , j - n + 1) / (j - n + 1.0) * gsl_pow_int(aux3, n);
        }
    } else if (t == 1) { /* a curve */
        tx = vf.x - vi.x;
        ty = vf.y - vi.y;
        l  = sqrt(tx * tx + ty * ty) ;
        tx = tx / l;
        ty = ty / l;
        theta = asin(0.5 * l / r);
        fori (nommt) {
            beta = theta * gp[2][i];
            xxgg = vf.x + r * ((cos(beta) - cos(theta)) *  ty + (sin(beta) - sin(theta)) * tx);
            yygg = vf.y + r * ((cos(beta) - cos(theta)) * -tx + (sin(beta) - sin(theta)) * ty);
            nxgg = cos(beta) * ty + sin(beta) * tx;
            wwgg = gw[2][i] * fabs(r * theta);
            aux1 = wwgg * nxgg;
            aux2 = xxgg - vr.x;
            aux3 = yygg - vr.y;
            c = 0 ;
            forj (o)
                forn (j +1)
                    m[c++] += aux1 * gsl_pow_int(aux2 , j - n + 1) / (j - n + 1.0) * gsl_pow_int(aux3, n);
        }
    }
}

flow compute_wall_condiction(border b, double inflow_angle, double xx, double yy, double Nxx, double Nyy) {
    if (b == INNER) {
        return WALL;       /* WALL CONDICTION */
    }
    else if (b == OUTTER) {
        if ((cos(inflow_angle)*Nxx + sin(inflow_angle)*Nyy) < 0)
            return INFLOW;   /* INFLOW condiction*/
        else
            return OUTFLOW;   /* OUTFLOW condiction */
    }
    return INTERNAL;
}

double * compute_radius_spline(int N ,double *tt ,double *xx, double *yy) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc () ;
    gsl_spline *spline_x = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline *spline_y = gsl_spline_alloc (gsl_interp_cspline, N);

    double *curve_radius = mdouble(N);

    gsl_spline_init(spline_x, tt, xx, N);
    gsl_spline_init(spline_y, tt, yy, N);

    int i ; 
    double dxi, dyi, d2xi, d2yi;

    for ( i = 0 ; i < N ; i++ ) {
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



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             MESHVIEWER INTEGRATION
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void v_draw_rawmesh(mesh * inm) {
    int i;

    if (inm->mesh_plot >= 0)
        return;

    float *v = (float *) malloc(sizeof(float) * inm->novertices * 2);
    fori (inm->novertices) {
        v[i * 2 + 0] = inm->vertices[i].x;
        v[i * 2 + 1] = inm->vertices[i].y;
    }

    unsigned int *id = (unsigned int *) malloc(sizeof(unsigned int) * inm->notriangles * 3);
    fori (inm->notriangles) {
        id[i * 3 + 0] = indexof(inm->vertices, inm->triangles[i].a);
        id[i * 3 + 1] = indexof(inm->vertices, inm->triangles[i].b);
        id[i * 3 + 2] = indexof(inm->vertices, inm->triangles[i].c);
    }

    //mv_add(MV_2D_TRIANGLES, v, inm->novertices, id, inm->notriangles * 3, mv_dgray, 0, NULL);

    int n;
    mv_add(MV_2D_TRIANGLES_AS_LINES, v, inm->novertices, id, inm->notriangles * 3, mv_white, 1.0, 1, &n);
    mv_setrotate(n, -inm->inflow_angle * 180.0 / PI);
    inm->mesh_plot = n;

    free(v);
    free(id);



    /*int faces = 2, fi = 0, ci = 0;
    v = (float *) malloc(sizeof(float) * faces * 2 * 2);
    float *c = (float *) malloc(sizeof(float) * faces * 3 * 2);

    v[fi++] = inm->edgeface[0][0].vi->x;    v[fi++] = inm->edgeface[0][0].vi->y;
    c[ci++] = 1.0; c[ci++] = 0.0; c[ci++] = 0.0;

    v[fi++] = inm->edgeface[0][0].vf->x;    v[fi++] = inm->edgeface[0][0].vf->y;
    c[ci++] = 1.0; c[ci++] = 1.0; c[ci++] = 0.0;




    v[fi++] = inm->edgeface[0][inm->noedgeface[0]-1].vi->x;    v[fi++] = inm->edgeface[0][inm->noedgeface[0]-1].vi->y;
    c[ci++] = 0.0; c[ci++] = 1.0; c[ci++] = 0.0;

    v[fi++] = inm->edgeface[0][inm->noedgeface[0]-1].vf->x;    v[fi++] = inm->edgeface[0][inm->noedgeface[0]-1].vf->y;
    c[ci++] = 0.0; c[ci++] = 1.0; c[ci++] = 1.0;


    fi = 0;
    printf("(%f, %f)---(%f, %f)\n", v[fi++], v[fi++], v[fi++], v[fi++]);
    printf("(%f, %f)---(%f, %f)\n", v[fi++], v[fi++], v[fi++], v[fi++]);

    mv_add(MV_2D_LINES | MV_USE_COLOUR_ARRAY, v, faces * 2, NULL, 0, c, 2.0, 2, &n);
    mv_setrotate(n, -inm->inflow_angle * 180.0 / PI);
    free(v);
    free(v);*/



    //v  = (float *) malloc(sizeof(float) * inm->noedgeface[i] * 2);
    /*id = (unsigned int *) malloc(sizeof(unsigned int) * inm->noedgeface[0] * 2);
    fori (inm->noedgeface[0]) {
        id[i * 2 + 0] = indexof(inm->vertices, inm->edgeface[0][i].vi);
        id[i * 2 + 1] = indexof(inm->vertices, inm->edgeface[0][i].vf);
    }

    mv_add(MV_2D_LINES, v, inm->novertices, id, inm->noedgeface[0] * 2, mv_red, 2.0f, NULL);*/
}


double tripol(double scale, double tri[][2]) {
    /* clipping */
    if (scale < 0.0) scale = 0.0;
    if (scale > 1.0) scale = 1.0;

    if (scale < tri[1][0])
        return ((tri[1][1] - tri[0][1]) / (tri[1][0] - tri[0][0])) * (scale - tri[0][0]) + tri[0][1];
    else if (scale >= tri[1][0] && scale <= tri[2][0]) return tri[1][1];
    else /* (scale > tri[2][0]) */
        return ((tri[3][1] - tri[2][1]) / (tri[3][0] - tri[2][0])) * (scale - tri[2][0]) + tri[2][1];
}

double transpol(double val, double m1, double m2, double min, double max) {
    double d1m1 = m1  - min, d2m1 = max - m1, d1m2 = m2  - min, d2m2 = max - m2;
    if (val < m1) return d1m2 / d1m1 * val;
    else          return d2m2 / d2m1 * (val-m1) + m2;
}

/* A bit more contrast */
double jetr[4][2] = { {0.475, 0.0}, {0.625, 1.0}, {0.875, 1.0}, {1.0,   0.6} };
double jetg[4][2] = { {0.1,   0.0}, {0.375, 1.0}, {0.625, 1.0}, {0.875, 0.0} };
double jetb[4][2] = { {0.0,   0.6}, {0.1,   1.0}, {0.35,  1.0}, {0.525, 0.0} };

/* MATLAB JET */
/*double jetr[4][2] = { {0.375, 0.0}, {0.625, 1.0}, {0.875, 1.0}, {1.0,   0.6} };
double jetg[4][2] = { {0.1,   0.0}, {0.375, 1.0}, {0.625, 1.0}, {0.875, 0.0} };
double jetb[4][2] = { {0.0,   0.6}, {0.1,   1.0}, {0.35,  1.0}, {0.625, 0.0} };*/


void v_draw_coefs(mesh * inm, int uselog) {
    #define max(a, b) ((a) > (b) ? (a) : (b))
    #define min(a, b) ((a) < (b) ? (a) : (b))

    int i;
    int id = 3; /* what to draw: 3 Cp, 0: Cm */
    double max = DBL_MIN, min = DBL_MAX, zero = 0.6;
    double dim = 1.0;

    fori (inm->novertices) {
        max = max(inm->u[i][id], max);
        min = min(inm->u[i][id], min);
        if (inm->vertices +i == inm->edgeface[1][inm->noedgeface[1] / 4].vi)
            zero = inm->u[i][id];
    }

    float *c = (float *) malloc(sizeof(float) * inm->novertices * 3);
    fori (inm->novertices) {
        double val = inm->u[i][id];
        val = (val - min) / (max - min);

        val = transpol(val, (zero - min) / (max - min), 0.5, 0.0, 1.0);

        if (uselog)
            val = log(1 + val) / log(2);

        c[i * 3 + 0] = tripol(val, jetr) * dim;
        c[i * 3 + 1] = tripol(val, jetg) * dim;
        c[i * 3 + 2] = tripol(val, jetb) * dim;
    }



    if ((inm->pressure_plot < 0 && !uselog) ||
        (inm->pressure_plotlog < 0 && uselog)) {

        float *v = (float *) malloc(sizeof(float) * inm->novertices * 2);
        fori (inm->novertices) {
            v[i * 2 + 0] = inm->vertices[i].x;
            v[i * 2 + 1] = inm->vertices[i].y;
        }

        unsigned int *id = (unsigned int *) malloc(sizeof(unsigned int) * inm->notriangles * 3);
        fori (inm->notriangles) {
            id[i * 3 + 0] = indexof(inm->vertices, inm->triangles[i].a);
            id[i * 3 + 1] = indexof(inm->vertices, inm->triangles[i].b);
            id[i * 3 + 2] = indexof(inm->vertices, inm->triangles[i].c);
        }

        int n;
        mv_add(MV_2D_TRIANGLES | MV_USE_COLOUR_ARRAY, v, inm->novertices, id, inm->notriangles * 3, c, 1.0, 0, &n);
        mv_setrotate(n, -inm->inflow_angle * 180.0 / PI);

        if (uselog)
            inm->pressure_plotlog = n;
        else
            inm->pressure_plot = n;

        if (inm->mesh_plot >= 0)
            mv_hide(inm->mesh_plot);

        free(v);
        free(id);
    } else {

        /*mv_setrotate(inm->pressure_plot, -inm->inflow_angle * 180.0 / PI);*/
        /*inm->inflow_angle -= 6 * PI / 180.0;*/

        if (uselog)
            mv_updatecolourarray(inm->pressure_plotlog, c, inm->novertices);
        else
            mv_updatecolourarray(inm->pressure_plot, c, inm->novertices);

    }

    /* VELOCITY */
    double factor = 0.0002;
    float *v = (float *) malloc(sizeof(float) * inm->novertices * 2 * 2);
    float *cv = (float *) malloc(sizeof(float) * inm->novertices * 3 * 2);
    fori (inm->novertices) {
        v[i * 4 + 0] = inm->vertices[i].x;
        v[i * 4 + 1] = inm->vertices[i].y;

        v[i * 4 + 2] = inm->vertices[i].x + inm->u[i][1] * factor;
        v[i * 4 + 3] = inm->vertices[i].y + inm->u[i][2] * factor;


        double r = c[i * 3 + 0];
        double g = c[i * 3 + 1];
        double b = c[i * 3 + 2];

        /*cv[i * 6 + 3] = cv[i * 6 + 0] = ((r + g + b) / 3.0) > 0.5 ? 0 : 1;
        cv[i * 6 + 4] = cv[i * 6 + 1] = ((r + g + b) / 3.0) > 0.5 ? 0 : 1;
        cv[i * 6 + 5] = cv[i * 6 + 2] = ((r + g + b) / 3.0) > 0.5 ? 0 : 1;*/
        cv[i * 6 + 3] = cv[i * 6 + 0] = 1 - r;
        cv[i * 6 + 4] = cv[i * 6 + 1] = 1 - g;
        cv[i * 6 + 5] = cv[i * 6 + 2] = 1 - b;
    }

    if (inm->velocity_plot < 0) {
        unsigned int *id = (unsigned int *) malloc(sizeof(unsigned int) * inm->novertices * 2);
        fori (inm->novertices * 2) {
            id[i] = i;
        }

        int n;
        mv_add(MV_2D_LINES | MV_USE_COLOUR_ARRAY, v, inm->novertices * 2, id, inm->novertices * 2, cv/*mv_pink*/, 1.0, 2, &n);
        mv_setrotate(n, -inm->inflow_angle * 180.0 / PI);
        inm->velocity_plot = n;

        free(id);
    }
    else {
        mv_updatevertexarray(inm->velocity_plot, v, inm->novertices * 2);
        mv_updatecolourarray(inm->velocity_plot, cv, inm->novertices * 2);
    }

    free(v);
    free(c);
    free(cv);
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             TIME MEASUREMENT
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
times * times_init(int initsize) {
    times * t = (times *) malloc(sizeof(times));

    t->tticks = (ttick *) malloc(sizeof(ttick) * initsize);
    t->size   = initsize;
    t->count   = 0;

    /* the initial tick */
    times_tick(t, "");
    t->tticks[t->count - 1].skip = -1;

    return t;
}

void times_tick(times * t, char * desc) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    /* if we reach the maximum size, allocate more! */
    if (t->count == t->size) {
        ttick    *ntt = (ttick *) malloc(sizeof(ttick) * t->size * 2);
        memcpy(ntt, t->tticks, sizeof(ttick) * t->size );
        free(t->tticks);
        t->tticks = ntt;

        t->size *= 2;
        //INFOM("Just to let you know that you've reached the maximum time ticks allocated. I've increased it from %d to %d", t->size / 2, t->size);
    }

    t->tticks[t->count].ticku = usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec / 10.0e6);
    t->tticks[t->count].ticks = usage.ru_stime.tv_sec + (usage.ru_stime.tv_usec / 10.0e6);

    if (desc == NULL)
        sprintf(t->tticks[t->count].desc, "step %d", t->count);
    /*else if (desc[0] == '\0')
        t->tticks[t->count].desc[0] = '\0';*/
    else
        sprintf(t->tticks[t->count].desc, "'%s'", desc);

    t->tticks[t->count].skip = 0;
    t->count++;
}

void times_zero(times * t) {
    times_tick(t, "");
    t->tticks[t->count - 1].skip = 1;
}

void times_print(times * t, int verbose) {
    int i;
    char tmp[256 + TIMEDESC];
    double elapsedu, elapseds,
           telapsedu = 0.0, telapseds = 0.0;

    /* compute the maximum string lengh */
    int maxmargin = 0;
    for (i = 1; i < t->count; i++) {
        int m = strlen(t->tticks[i].desc);
            maxmargin = m > maxmargin ? m : maxmargin;
    }
    maxmargin += 8;



    /* count total */
    for (i = 1; i < t->count; i++) {
        if (t->tticks[i].skip) {
            continue;
        }

        elapsedu = t->tticks[i].ticku - t->tticks[i - 1].ticku;
        elapseds = t->tticks[i].ticks - t->tticks[i - 1].ticks;

        telapsedu += elapsedu;
        telapseds += elapseds;
    }

    int maxd = 6;
    int maxs = (int) log10(telapsedu+telapseds) + 1 + maxd + 1;

    if (verbose)
        for (i = 1; i < t->count; i++) {
            if (t->tticks[i].skip) {
                continue;
            }

            elapsedu = t->tticks[i].ticku - t->tticks[i - 1].ticku;
            elapseds = t->tticks[i].ticks - t->tticks[i - 1].ticks;

            sprintf(tmp, "Time in %s", t->tticks[i].desc);
            //printf("%-*s: %fs (user %fs + sys %fs)\n", maxmargin, tmp, elapsedu + elapseds, elapsedu, elapseds);
            INFOM("%-*s: %*.*fs (user %*.*fs + sys %*.*fs)", maxmargin, tmp, maxs, maxd, elapsedu + elapseds, maxs, maxd, elapsedu, maxs, maxd, elapseds);
        }

    /* Total time */
    sprintf(tmp, "Total time");
    //printf("%-*s: %fs (user %fs + sys %fs)\n", maxmargin, tmp, elapsedu + elapseds, elapsedu, elapseds);
    INFOM(BOLD"%-*s"RESET": "BOLD"%*.*fs"RESET" (user %*.*fs + sys %*.*fs)", maxmargin, tmp, maxs, maxd, telapsedu + telapseds, maxs, maxd, telapsedu, maxs, maxd, telapseds);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             Signal Handling
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void sig_handler(int signo) {
    /*
     * This handle signals
     */
    //static char *FUN = "";
    if (signo == SIGINT) {
        printf("\r");
        fprintf(stderr, "\r");
        /* VISUALIZER */
        //mv_stop();

        if (!running_solver) {
            ERRM("SIGINT received. Aborting!");
            exit(0);
        }

        if (interrupt_solver && running_solver) {
            ERRM("SIGINT received. Forcing abortion!");
            exit(0);
        }

        interrupt_solver = 1;
        INFOM("SIGINT received. Interrupting the solver as soon as we can!");
    }
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                             DEBUG
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void print2d(char * p, double **u, int l, int c, int new) {
    FILE *f;

    if (new) f = fopen("debug", "w");
    else     f = fopen("debug", "a");

    int i, j;
    for (i = 0; i < l; i++) {
        fprintf(f, "[%s] %4d: ", p , i);
        for (j = 0; j < c; j++)
            fprintf (f, "%.50e ", u[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);
}
void print3d(double ***u, int m, int l, int c, int new) {
    FILE *f;

    if (new) f = fopen("debug", "w");
    else     f = fopen("debug", "a");

    int i, j, n;
    forn (m)
    for (i = 0; i < l; i++) {
        fprintf(f, "[coef] %4d: ", i);
        for (j = 0; j < c; j++)
            fprintf (f, "%.50e ", u[n][i][j]);
        fprintf(f, "\n");
    }
        fclose(f);
}





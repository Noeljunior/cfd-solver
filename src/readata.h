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


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                          DATA STRUCTURE
 *                                           and typedefs
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  The struct to return the vertices, edges, and triangles
 */
typedef struct cfdrd_ds {
    double      **vertices;         /* sizev * : vertices [{ xi, yi, bi, ri }] */
//    int         *border;            /* sizev * : border [ bi ] */
//    double      *radius;            /* sizev * : vertices [ ri ] */
    int         sizev;              /* vertices count */
    //int         trailingid;         /* the vertice id of the trailing vertice */

    int         **triangles;        /* sizet * 3 : edges [{ v1, v2, v3 }] */
    int         sizet;


    double      **params;           /* sizep * : parametric function's parameters
                                        [{ border_id, func_id, #params, param * #params }] */
    int         sizep;              /* count of params functions */

    double      chord;              /* length of the airfoil's chord */

    int         **edges;            /* sizee * 3 : edges [{ v1, v2, border_type }] */
    int         sizee;
} cfdrd_ds;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *                                       FUNCTION DECLARATION
 *                                       of external interface
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Try to find out the version of file and read it accordingly
 */
cfdrd_ds *      cfdrd_readfile_auto(char *path, int rstdin, int force, int bequiet);
void            cfdrd_free(cfdrd_ds * ds);



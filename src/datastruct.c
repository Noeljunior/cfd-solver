#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#include "datastruct.h"

/* Print info related to this file */
#define INFO(i) printf("INFO[DS]: %s\n", i)



int ds_detect_file_type(FILE * input) {
    /* TODO: read first line(s)
             detect which file type are we talking about
             seek to the begining of the file
             return that type
     */
}


/*
 * Read the input file content
 */
mesh * ds_read_file(FILE *fid) {
    mesh * in = (mesh *) malloc(sizeof(mesh));
    int i;
    
    assert(fid == NULL, "Could not open input file.");
    
    #ifdef OLD                      /* The OLD input file */
    char trash[1024];
    int  trashi;
    fscanf(fid, "%s %d", trash, &trashi);           /* trash */
    fscanf(fid, "%s %d", trash, &in->dimension);     /* dimension */
    fscanf(fid, "%s %d", trash, &in->novertices);    /* number of vertices */
    in->vertices = (coord *) malloc(sizeof(coord) * in->novertices);
    for (i = 0; i < in->novertices; i++)            /* vertices */
        fscanf(fid, "%lf %lf %d %d", &in->vertices[i].x, &in->vertices[i].y, &trashi, &trashi);
    
    fscanf(fid, "%s %d", trash, &in->noedges);       /* number of edges */
    in->edges = (edge *) malloc(sizeof(edge) * in->noedges);
    for (i = 0; i < in->noedges; i++) {              /* edges */
        fscanf(fid, "%d %d %d", &in->edges[i].a, &in->edges[i].b, &in->edges[i].border);
        in->edges[i].a--;
        in->edges[i].b--;
        in->edges[i].border--;
    }
    
    fscanf(fid, "%s %d", trash, &in->notriangles);   /* number of triangles */
    in->triangles = (triangle *) malloc(sizeof(triangle) * in->notriangles);
    for (i = 0; i < in->notriangles; i++) {         /* triangles */
        fscanf(fid, "%d %d %d %d", &in->triangles[i].a, &in->triangles[i].b, &in->triangles[i].c, &trashi);
        in->triangles[i].a--;
        in->triangles[i].b--;
        in->triangles[i].c--;
    }
    
    fscanf(fid, "%s", trash);                       /* EOF */
    
    
    #else                           /* The new input file */
    in->dimension = 3;
    fscanf(fid, "%d", &in->novertices);              /* number of vertices */
    in->vertices = (coord *) malloc(sizeof(coord) * in->novertices);
    for (i = 0; i < in->novertices; i++)            /* vertices */
        fscanf(fid, "%lf %lf", &in->vertices[i].x, &in->vertices[i].y);
    
    fscanf(fid, "%d", &in->noedges);                 /* number of edges */
    in->edges = (edge *) malloc(sizeof(edge) * in->noedges);
    for (i = 0; i < in->noedges; i++)               /* edges */
        fscanf(fid, "%d %d %d", &in->edges[i].a, &in->edges[i].b, &in->edges[i].border);
    
    fscanf(fid, "%d", &in->notriangles);             /* number of triangles */
    in->triangles = (triangle *) malloc(sizeof(triangle) * in->notriangles);
    for (i = 0; i < in->notriangles; i++)           /* triangles */
        fscanf(fid, "%d %d %d", &in->triangles[i].a, &in->triangles[i].b, &in->triangles[i].c);
    #endif
    
    return in;
}





/* Abort the execution saying why by text and errno */
void abortc(char *why, int errno) {
    fprintf(stderr, "Err: %s\n", why); exit(errno);
}
/* Abort the execution saying why by text */
void abortw(char *why) {
    abortc(why, -1);
}
/* Abort the execution */
void aborte() {
    abortc("Something went wrong!", -1);
}
/* Abort if condition is true */
void assert(int cond, char * why) {
    if (cond) abortc(why, -1);
}
/* Abort if condition is false */
void nassert(int cond, char * why) {
    if (!cond) abortc(why, -1);
}

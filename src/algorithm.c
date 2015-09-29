#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#include "datastruct.h"
#include "algorithm.h"

/* Print info related to this file */
#define INFO(i) printf("INFO[AL]: %s\n", i)

#define fori(cond) for (i = 0; i < cond; i++)
#define forj(cond) for (j = 0; j < cond; j++)

void ** malloc2d(size_t s, int n, int m) {
    int i;
    void ** res = malloc(sizeof(void *) * n);
    
    void * t = malloc(s * n * m);
    
    for (i = 0; i < n; i++) {
        res[i] = t + i * m;
    }
    
    return res;
}

void compute(args * ina, mesh * inm) {
    int i;

    inm->novt   = inm->novertices + inm->notriangles;
    inm->noc    = (ina->ORDER * (ina->ORDER + 1)) / 2;
    inm->nommt  = (int) (ceil((ina->ORDER + 1.0) / 2.0));
    inm->lpd    = (int) (ceil(ina->ORDER / 2.0 ));
    inm->Aexp   = 1;
    inm->norv   = 4;
    inm->kdelta = 0.25;
    inm->m1     = 0.8;
    inm->m2     = 0.85;

    inm->faces = (face *) malloc(sizeof(face) * inm->novt);

    for (i = 0; i < 3; i++) {
        inm->gps[i] = (double *) malloc(sizeof(double) * (i+1));
        inm->gws[i] = (double *) malloc(sizeof(double) * (i+1));
    }
    
    inm->tcentroids = (coord *) malloc(sizeof(coord) * inm->notriangles);
    
    inm->mmt = (double **) malloc2d(sizeof(double), inm->novertices, inm->noc);
    
    //TODO continue on AVC LeituraArqyuivo.c:114
}

face * compute_faces(args * ina, mesh * inm) {
    face * faces = (face *) malloc(sizeof(face) * inm->novt);
    
    int i, j, k, p, t;
    
    int pos[inm->notriangles];
    
    fori (inm->notriangles) pos[i] = i; // TODO i or i + 1
    
    
    /* TODO
     * these two for's, for each edge, find a triangle that have that edge
     * in other words, this could be done by an hash table
     */
    /*p = 0;
    fori (inm->noedges) {
        
        forj (inm->notriangles) {
            k = pos[j];
            
            if (inm->edges[i].a == inm->triangles[j].a &&
                inm->edges[i].b == inm->triangles[j].b) {
                face[p].vi = inm->edges[i].a;
                face[p].vf = inm->edges[i].b;
                face[p].tl = k;
                face[p].tr = 0;
                face[p].vl = inm->triangles[j].c;
                face[p].vr = 0;
                p++;
                break;
            }
            if (inm->edges[i].a == inm->triangles[j].b &&
                inm->edges[i].b == inm->triangles[j].c) {
                face[p].vi = inm->edges[i].a;
                face[p].vf = inm->edges[i].b;
                face[p].tl = k;
                face[p].tr = 0;
                face[p].vl = inm->triangles[j].a;
                face[p].vr = 0;
                p++;
                break;
            }
            if (inm->edges[i].a == inm->triangles[j].c &&
                inm->edges[i].b == inm->triangles[j].a) {
                face[p].vi = inm->edges[i].a;
                face[p].vf = inm->edges[i].b;
                face[p].tl = k;
                face[p].tr = 0;
                face[p].vl = inm->triangles[j].b;
                face[p].vr = 0;
                p++;
                break;
            }
        }
        
        t      = pos[j];
        pos[j] = pos[i];
        pos[i] = t;
        
    }
    */
    
    
    return faces;
}
























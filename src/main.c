#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cfdsolver.h"
#include "readata.h"
#include "ui.h"

#include "infocolour.h"

int main(int argc, char ** argv) {

    ui_args *args = ui_readargs(argc, argv);


    cfds_args * inargs = (cfds_args *) malloc(sizeof(cfds_args));
    inargs->order          = args->order;
    inargs->cfl            = args->cfl;
    inargs->max_iterations = args->max_iterations;
    inargs->nr_threashold  = args->nr_threashold;
    inargs->verbose        = args->verbose;
    inargs->quiet          = args->quiet;
    inargs->showinner      = args->showinner;
    inargs->fclassify      = args->fclassify;


    cfdrd_ds * ds = cfdrd_readfile_auto(args->mesh_file, args->rstdin, args->mftype, args->quiet);

    cfds_mesh *m = cfds_init(inargs, ds->vertices, ds->sizev, ds->edges, ds->sizee, ds->triangles, ds->sizet);

    cfdrd_free(ds);

    cfds_solve(m);
    cfds_free(m);




}









#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <stddef.h>

#include "common.h"
#include "body.h"
#include "cell.h"
#include "darray.h"

#define MASTER 0
#define MAX_LINE_LENGTH 1024

void get_universe_size(vect3_t universe_min, vect3_t universe_max, darray_t* bodies) {
    universe_min[0] = universe_min[1] = universe_min[2] = DBL_MAX;
    universe_max[0] = universe_max[1] = universe_max[2] = DBL_MIN;

    for (int i = 0; i < bodies->length; i++) {
        body_t* body = AS_BODY_PTR(darray_get(bodies, i));
        for (int j = 0; j < DIMS; j++) {
            if (body->position[j] < universe_min[j]) {
                universe_min[j] = body->position[j];
            }
            if (body->position[j] > universe_max[j]) {
                universe_max[j] = body->position[j];
            }
        }
    }

    for(int c = 0; c < DIMS; c++){
        double global_min, global_max;
        MPI_Allreduce(&universe_min[c], &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&universe_max[c], &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        universe_min[c] = global_min;
        universe_max[c] = global_max;
    }
}

int main(int argc, char const *argv[]) {
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int n_bodies = atoi(argv[2]);
    int iterations = atoi(argv[3]);

    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    init_MPI_Body_datatype();
    init_MPI_Cell_datatype();

    // Scatter bodies across ranks
    int* sendcounts = malloc(world_size * sizeof(int));
    int* displs = malloc(world_size * sizeof(int));
    
    int block_size = (int) (n_bodies + world_size - 1) / world_size;
    for (int i = 0; i < world_size; i++) {
        int start = i * block_size;
        int end = MIN(start + block_size, n_bodies);        
        sendcounts[i] = end - start;
        displs[i] = start;
    }

    body_t* bodies = NULL;

    if (rank == MASTER) {
        // Initialize and read the data structures
        bodies = (body_t*) malloc(n_bodies * sizeof(body_t));

        char str[MAX_LINE_LENGTH];
        for (int i = 0; i < n_bodies; i++) {
            double x, y, z, vx, vy, vz, mass;
            if (fgets(str, MAX_LINE_LENGTH, f) == NULL) {
                printf("Error reading file.\n");
                exit(1);
            }
            sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz, &mass);
            bodies[i].position[X] = x;
            bodies[i].position[Y] = y;
            bodies[i].position[Z] = z;

            bodies[i].velocity[X] = vx;
            bodies[i].velocity[Y] = vy;
            bodies[i].velocity[Z] = vz;

            bodies[i].mass = mass;
            bodies[i].id = i;
        }
    }
    int my_bodies_count = sendcounts[rank];
    body_t* my_bodies = malloc(my_bodies_count * sizeof(body_t));
    MPI_Scatterv(
        bodies,
        sendcounts,
        displs,
        MPI_Body,

        my_bodies,
        my_bodies_count,
        MPI_Body,
        MASTER,
        MPI_COMM_WORLD
    );
    // Done initial scatter

    for (int i = 0; i < my_bodies_count; i++) {
        printf("rank %d id %d %f\n", rank, my_bodies[i].id, my_bodies[i].position[X]);
    }

    darray_t* dbodies = malloc(sizeof(darray_t));
    darray_init(dbodies, 1024);

    for (int i = 0; i < my_bodies_count; i++) {
        printf("rank %d appending %d\n", rank, my_bodies[i].id); fflush(stdout);
        darray_append(dbodies, AS_VOID_PTR(&my_bodies[i]));
    }

    printf("dbodies len: %zu\n", dbodies->length);
    for (int i = 0; i < dbodies->length; i++) {
        body_t* body = AS_BODY_PTR(dbodies->elements[i]);
        printf("rank %d id %d %f\n", rank, body->id, body->position[X]);
    }

    vect3_t universe_min, universe_max;
    get_universe_size(universe_min, universe_max, dbodies);

    if (rank == MASTER) {
        printf("%f %f\n", universe_min[X], universe_max[X]);
        printf("%f %f\n", universe_min[Y], universe_max[Y]);
        printf("%f %f\n", universe_min[Z], universe_max[Z]);
    }

    

    MPI_Finalize();
}

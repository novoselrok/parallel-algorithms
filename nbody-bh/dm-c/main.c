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
#include "orb.h"

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

int main(int argc, char *argv[]) {
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
    
    int block_size = (n_bodies + world_size - 1) / world_size;
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
            bodies[i].work = 1;
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

    darray_t* dbodies = malloc(sizeof(darray_t));
    darray_init(dbodies, 1024);

    for (int i = 0; i < my_bodies_count; i++) {
        darray_append(dbodies, AS_VOID_PTR(&my_bodies[i]));
    }
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    double dt = 0.1;
    size_t n_splits = (size_t) log2(world_size);
    darray_t* my_bounds = malloc(sizeof(darray_t));
    darray_t* other_bounds = malloc(sizeof(darray_t));
    darray_t* partners = malloc(sizeof(darray_t));
    darray_init(my_bounds, n_splits);
    darray_init(other_bounds, n_splits);
    darray_init(partners, n_splits);

    for (int iter = 0; iter < iterations; iter++) {
        vect3_t universe_min, universe_max;
        get_universe_size(universe_min, universe_max, dbodies);

        darray_clear(my_bounds);
        darray_clear(other_bounds);
        darray_clear(partners);

        orb(dbodies, rank, world_size, universe_min, universe_max, my_bounds, other_bounds, partners);

        cell_t* root = malloc(sizeof(cell_t));
        init_cell_with_bounds(root, universe_min, universe_max);

        for (size_t i = 0; i < my_bounds->length; i++) {
            double* tuple = AS_DOUBLE_PTR(darray_get(my_bounds, i));
            vect3_t min, max;
            unpack_tuple(tuple, min, max);
            insert_empty_cell(root, min, max);
        }

        for (size_t i = 0; i < dbodies->length; i++) {
            body_t* body = AS_BODY_PTR(darray_get(dbodies, i));
            insert_body(root, body);
        }

        //printf("rank %d root.cm = %f %f %f\n", rank, root->cm[X], root->cm[Y], root->cm[Z]);

        darray_t* cells_to_send = malloc(sizeof(darray_t));
        darray_init(cells_to_send, 1024);
        MPI_Status status;
        for (size_t i = 0; i < my_bounds->length; i++) {
            vect3_t my_min, my_max, other_min, other_max;
            double* my_bounds_tuple = AS_DOUBLE_PTR(darray_get(my_bounds, i));
            double* other_bounds_tuple = AS_DOUBLE_PTR(darray_get(other_bounds, i));
            unpack_tuple(my_bounds_tuple, my_min, my_max);
            unpack_tuple(other_bounds_tuple, other_min, other_max);

            darray_clear(cells_to_send);
            get_cells_to_send(root, NULL, other_min, other_max, (int) i, 0, cells_to_send, rank);
            //printf("rank %d cells_to_send_len %d\n", rank, cells_to_send->length);
            mpi_cell_t* packed = pack_cells(cells_to_send);

            int* partner_above = AS_INT_PTR(darray_get(partners, i));
            int partner = partner_above[0];
            bool above_split = (bool) partner_above[1];

            int recv_cells_count;
            mpi_cell_t* recv_cells = NULL;
            if (above_split) {
                MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_Cell, &recv_cells_count);
                recv_cells = malloc(sizeof(mpi_cell_t) * recv_cells_count);
                MPI_Recv(recv_cells, recv_cells_count, MPI_Cell, partner, 0, MPI_COMM_WORLD, &status);

                MPI_Send(packed, (int) cells_to_send->length, MPI_Cell, partner, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(packed, (int) cells_to_send->length, MPI_Cell, partner, 0, MPI_COMM_WORLD);

                MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_Cell, &recv_cells_count);
                recv_cells = malloc(sizeof(mpi_cell_t) * recv_cells_count);
                MPI_Recv(recv_cells, recv_cells_count, MPI_Cell, partner, 0, MPI_COMM_WORLD, &status);
            }
            free(packed);

            darray_t* root_cells = reconstruct_received_cells(recv_cells, recv_cells_count);

            //printf("rank %d root_cells_len %d\n", rank, root_cells->length);

            for (size_t j = 0; j < root_cells->length; j++) {
                cell_t* root_cell = AS_CELL_PTR(darray_get(root_cells, j));
                insert_cell(root, root_cell);
            }

            free(recv_cells);
            darray_free(root_cells);
        }
        darray_free(cells_to_send);
        //printf("rank %d root.cm = %f %f %f root.mass = %f\n", rank, root->cm[X], root->cm[Y], root->cm[Z], root->mass);

        for (size_t i = 0; i < dbodies->length; i++) {
            body_t* body = AS_BODY_PTR(darray_get(dbodies, i));
            double start_force_time = MPI_Wtime();
            compute_force(root, body);
            body->work = MPI_Wtime() - start_force_time;
        }

        for (size_t i = 0; i < dbodies->length; i++) {
            body_t* body = AS_BODY_PTR(darray_get(dbodies, i));
            body->position[X] += dt * body->velocity[X];
            body->position[Y] += dt * body->velocity[Y];
            body->position[Z] += dt * body->velocity[Z];

            body->velocity[X] += dt / body->mass * body->force[X];
            body->velocity[Y] += dt / body->mass * body->force[Y];
            body->velocity[Z] += dt / body->mass * body->force[Z];
        }

        free_cell(root);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == MASTER) {
        printf("%f\n", end_time - start_time);
    }

    char name[50];
    sprintf(name, "out%d.txt", rank);
    FILE* outf;
    outf = fopen(name, "w");
    if(outf == NULL) {
        printf("Could not open output file.\n");
        exit(1);
    }
    for (size_t i = 0; i < dbodies->length; i++) {
        body_t* body = AS_BODY_PTR(darray_get(dbodies, i));
        fprintf(outf, "%d %e %e %e %e %e %e\n", body->id, body->position[X], body->position[Y], body->position[Z], body->velocity[X], body->velocity[Y], body->velocity[Z]);
    }
    fclose(outf);

    MPI_Finalize();
}

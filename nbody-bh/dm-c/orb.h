#ifndef _ORB_H
#define _ORB_H

#include <mpi.h>
#include <string.h>
#include <stdbool.h>

#include "common.h"
#include "darray.h"
#include "body.h"

#define BISECTION_MAX_ITER 200
#define BISECTION_TOL 1e-10

bool is_above_split(int rank, int n_procs_left) {
    int rel_bit = (int) log2(n_procs_left);
    int is_above = (rank >> (rel_bit - 1)) & 1; // get last bit
    return is_above == 1;
}

int get_partner_rank(int rank, int n_procs_left) {
    int rel_bit = (int) log2(n_procs_left);
    return rank ^ (1 << (rel_bit - 1)); // rank xor (2^(rel_bit-1))
}

double frac_weight_below(darray_t* bodies, double split, coord_t coord, MPI_Comm comm) {
    double work_below = 0.0, work_above = 0.0;
    for (size_t i = 0; i < bodies->length; i++) {
        body_t* body = AS_BODY_PTR(darray_get(bodies, i));

        if(body->position[coord] > split){
            work_above += body->work;
        } else {
            work_below += body->work;
        }
    }

    double all_work_below = 0.0, all_work_above = 0.0;
    MPI_Allreduce(&work_above, &all_work_above, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&work_below, &all_work_below, 1, MPI_DOUBLE, MPI_SUM, comm);

    double fraction = all_work_below / (all_work_above + all_work_below);

    return fraction - 0.5; // Subtract 0.5 to look for the median
}

double bisection(double min, double max, darray_t* bodies, coord_t coord, MPI_Comm comm) {
    double fmin = frac_weight_below(bodies, min, coord, comm);
    double fmid, mid = 0.0;

    int iter = 0;
    while (iter < BISECTION_MAX_ITER && fabs((max - min) / 2) > BISECTION_TOL) {
        mid = (min + max) / 2;
        fmid = frac_weight_below(bodies, mid, coord, comm);

        if (fmid == 0.0) {
            break;
        } else if (fmin * fmid > 0) {
            min = mid;
        } else {
            max = mid;
        }

        iter++;
    }

    return mid;
}

void orb(
        darray_t* bodies,
        int rank,
        int world_size,
        vect3_t universe_min,
        vect3_t universe_max,
        darray_t* my_bounds,
        darray_t* other_bounds,
        darray_t* partners) {

    MPI_Status status;
    MPI_Comm comm_subset;
    vect3_t my_min, my_max;
    memcpy(my_min, universe_min, sizeof(vect3_t));
    memcpy(my_max, universe_max, sizeof(vect3_t));

    darray_t* my_bodies = malloc(sizeof(darray_t));
    darray_t* other_bodies = malloc(sizeof(darray_t));
    darray_init(my_bodies, bodies->length);
    darray_init(other_bodies, bodies->length);

    int group = 0;
    bool above_split = false;

    size_t n_splits = (size_t) log2(world_size);
    for (int i = 0; i < n_splits; i++) {
        int n_procs_left = world_size / (int) pow(2, i);

        group = group << 1; // group * 2
        if (above_split) {
            group = group | 1; // set least significant bit to 1
        }

        // Split communication according to the group of each process
        MPI_Comm_split(MPI_COMM_WORLD, group, rank, &comm_subset);

        // At each split select the next coordinate on which to split
        coord_t coord = i % DIMS;

        double split = bisection(my_min[coord], my_max[coord], bodies, coord, comm_subset);

        //printf("rank %d iter %d split %f\n", rank, i, split);

        MPI_Comm_free(&comm_subset);

        above_split = is_above_split(rank, n_procs_left);

        vect3_t other_min, other_max;
        memcpy(other_min, my_min, sizeof(vect3_t));
        memcpy(other_max, my_max, sizeof(vect3_t));

        // Resize my bounds and other bounds
        if (above_split) {
            my_min[coord] = split;
            other_max[coord] = split;
        } else {
            my_max[coord] = split;
            other_min[coord] = split;
        }

        //printf("rank %d %f %f %f %f %f %f\n", rank, other_min[X], other_min[Y], other_min[Z], other_max[X], other_max[Y], other_max[Z]);

        // Save my bounds and other bounds
        double* my_bounds_tuple = malloc(sizeof(double) * DIMS * 2);
        double* other_bounds_tuple = malloc(sizeof(double) * DIMS * 2);
        pack_tuple(my_bounds_tuple, my_min, my_max);
        pack_tuple(other_bounds_tuple, other_min, other_max);

        darray_append(my_bounds, AS_VOID_PTR(my_bounds_tuple));
        darray_append(other_bounds, AS_VOID_PTR(other_bounds_tuple));

        darray_clear(my_bodies);
        darray_clear(other_bodies);
        for (size_t body_idx = 0; body_idx < bodies->length; body_idx++) {
            body_t* body = AS_BODY_PTR(darray_get(bodies, body_idx));
            if ((body->position[coord] - split > 0) == above_split) {
                darray_append(my_bodies, AS_VOID_PTR(body));
            } else {
                darray_append(other_bodies, AS_VOID_PTR(body));
            }
        }

        //printf("rank %d iter %d my_bodies_len %d\n", rank, i, my_bodies->length);

        int partner = get_partner_rank(rank, n_procs_left);
        int* partner_above = (int*) malloc(sizeof(int) * 2);
        partner_above[0] = partner;
        partner_above[1] = above_split;
        darray_append(partners, AS_VOID_PTR(partner_above));

        body_t* recv_bodies = NULL;
        int recv_bodies_count;
        body_t* packed = pack_bodies(other_bodies);
        if (above_split) {
            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_Body, &recv_bodies_count);
            recv_bodies = malloc(sizeof(body_t) * recv_bodies_count);
            MPI_Recv(recv_bodies, recv_bodies_count, MPI_Body, partner, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            MPI_Send(packed, (int) other_bodies->length, MPI_Body, partner, 0, MPI_COMM_WORLD);
        } else {
            MPI_Send(packed, (int) other_bodies->length, MPI_Body, partner, 0, MPI_COMM_WORLD);

            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_Body, &recv_bodies_count);
            recv_bodies = malloc(sizeof(body_t) * recv_bodies_count);
            MPI_Recv(recv_bodies, recv_bodies_count, MPI_Body, partner, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
        free(packed);

        for (size_t body_idx = 0; body_idx < recv_bodies_count; body_idx++) {
            body_t* body = &recv_bodies[body_idx];
            darray_append(my_bodies, AS_VOID_PTR(body));
        }

        darray_copy(bodies, my_bodies);
    }

    //printf("rank %d final_bodies_len %d\n", rank, bodies->length);
}

#endif

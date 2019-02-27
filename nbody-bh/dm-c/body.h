#ifndef _BODY_H
#define _BODY_H

#include <mpi.h>
#include <stddef.h>

#include "common.h"

#define AS_BODY_PTR(x) ((body_t*)(x))

typedef struct _body_t {
    int id;
    vect3_t force;
    vect3_t position;
    vect3_t velocity;
    double mass;
} body_t;

extern MPI_Datatype MPI_Body;

MPI_Datatype MPI_Body;

void init_MPI_Body_datatype() {
    int blocklengths[5] = {1, 3, 3, 3, 1};
    MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[5];

    offsets[0] = offsetof(body_t, id);
    offsets[1] = offsetof(body_t, force);
    offsets[2] = offsetof(body_t, position);
    offsets[3] = offsetof(body_t, velocity);
    offsets[4] = offsetof(body_t, mass);

    MPI_Type_create_struct(5, blocklengths, offsets, types, &MPI_Body);
    MPI_Type_commit(&MPI_Body);
}

void reset_force(body_t* body) {
    body->force[X] = 0.0;
    body->force[Y] = 0.0;
    body->force[Z] = 0.0;
}

#endif

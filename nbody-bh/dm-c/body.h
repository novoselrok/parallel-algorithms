#ifndef _BODY_H
#define _BODY_H

#include <mpi.h>
#include <stddef.h>

#include "common.h"
#include "darray.h"

#define AS_BODY_PTR(x) ((body_t*)(x))

typedef struct _body_t {
    int id;
    vect3_t force;
    vect3_t position;
    vect3_t velocity;
    double mass;
    double work;
} body_t;

extern MPI_Datatype MPI_Body;

MPI_Datatype MPI_Body;

void init_MPI_Body_datatype() {
    int blocklengths[6] = {1, 3, 3, 3, 1, 1};
    MPI_Datatype types[6] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[6];

    offsets[0] = offsetof(body_t, id);
    offsets[1] = offsetof(body_t, force);
    offsets[2] = offsetof(body_t, position);
    offsets[3] = offsetof(body_t, velocity);
    offsets[4] = offsetof(body_t, mass);
    offsets[5] = offsetof(body_t, work);

    MPI_Type_create_struct(6, blocklengths, offsets, types, &MPI_Body);
    MPI_Type_commit(&MPI_Body);
}

body_t* pack_bodies(darray_t* darray) {
    body_t* packed = malloc(darray->length * sizeof(body_t));
    for (size_t i = 0; i < darray->length; i++) {
        body_t* body = AS_BODY_PTR(darray_get(darray, i));
        memcpy(&packed[i], body, sizeof(body_t));
    }
    return packed;
}

void reset_force(body_t* body) {
    body->force[X] = 0.0;
    body->force[Y] = 0.0;
    body->force[Z] = 0.0;
}

#endif

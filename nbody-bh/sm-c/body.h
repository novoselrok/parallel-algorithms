#ifndef _BODY_H
#define _BODY_H

#include "common.h"

typedef struct _body_t {
    int id;
    vect3_t force;
    vect3_t position;
    vect3_t velocity;
    double mass;
} body_t;

void reset_force(body_t* body) {
    body->force[X] = 0.0;
    body->force[Y] = 0.0;
    body->force[Z] = 0.0;
}

#endif

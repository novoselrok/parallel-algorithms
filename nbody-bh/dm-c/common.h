#ifndef _COMMON_H
#define _COMMON_H

#include <math.h>
#include <string.h>

#define G 6.67e-11
#define DIMS 3
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define AS_DOUBLE_PTR(x) ((double*)(x))
#define AS_INT_PTR(x) ((int*)(x))

typedef enum _coord_t {
    X = 0,
    Y = 1,
    Z = 2
} coord_t;

typedef double vect3_t[3];

void init_vect(vect3_t vec) {
    vec[X] = 0.0;
    vec[Y] = 0.0;
    vec[Z] = 0.0;
}

void pack_tuple(double* tuple, vect3_t min, vect3_t max) {
    memcpy(tuple, min, sizeof(vect3_t));
    memcpy(&tuple[3], max, sizeof(vect3_t));
}

void unpack_tuple(double* tuple, vect3_t min, vect3_t max) {
    memcpy(min, tuple, sizeof(vect3_t));
    memcpy(max, &tuple[3], sizeof(vect3_t));
}

double distance(vect3_t vec1, vect3_t vec2) {
    double dx = vec1[X] - vec2[X];
    double dy = vec1[Y] - vec2[Y];
    double dz = vec1[Z] - vec2[Z];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

#endif

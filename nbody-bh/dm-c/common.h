#ifndef _COMMON_H
#define _COMMON_H

#include <math.h>

#define G 6.67e-11
#define DIMS 3
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

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

double distance(vect3_t vec1, vect3_t vec2) {
    double dx = vec1[X] - vec2[X];
    double dy = vec1[Y] - vec2[Y];
    double dz = vec1[Z] - vec2[Z];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

#endif

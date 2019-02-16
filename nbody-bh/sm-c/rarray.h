#ifndef _RARRAY_H
#define _RARRAY_H

#include <stdlib.h>
#include <stdio.h>

#include "body.h"

typedef struct _rarray_t {
    body_t** bodies;
    size_t length;
    size_t capacity;
} rarray_t;


void rarray_init(rarray_t* rarray, size_t capacity) {
    rarray->bodies = (body_t**) malloc(capacity * sizeof(body_t*));
    rarray->capacity = capacity;
    rarray->length = 0;
}

void rarray_resize(rarray_t* rarray, size_t new_capacity) {
    printf("Resizing to %d\n", new_capacity); fflush(stdout);
    rarray->bodies = realloc(rarray->bodies, new_capacity * sizeof(body_t*));
}

void rarray_append(rarray_t* rarray, body_t* body) {
    if (rarray->length == rarray->capacity) {
        rarray->capacity *= 2;
        rarray_resize(rarray, rarray->capacity);
    }

    rarray->bodies[rarray->length] = body;
    rarray->length++;
}

// TODO: Convert to macro
body_t* rarray_get(rarray_t* rarray, size_t idx) {
    return rarray->bodies[idx];
}

void rarray_free(rarray_t* rarray) {
    for (int i = 0; i < rarray->capacity; i++) {
        free(rarray->bodies[i]);
    }
    free(rarray->bodies);
    free(rarray);
}

#endif

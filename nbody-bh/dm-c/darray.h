#ifndef _DARRAY_H
#define _DARRAY_H

#include <stdlib.h>
#include <stdio.h>

typedef struct _darray_t {
    void** elements;
    size_t length;
    size_t capacity;
} darray_t;

#define AS_VOID_PTR(x) ((void*)(x))

void darray_init(darray_t* darray, size_t capacity) {
    darray->elements = (void**) malloc(capacity * sizeof(void*));
    darray->capacity = capacity;
    darray->length = 0;
}

void darray_resize(darray_t* darray, size_t new_capacity) {
    darray->elements = (void**) realloc(darray->elements, new_capacity * sizeof(void*));
}

void darray_append(darray_t* darray, void* element) {
    if (darray->length == darray->capacity) {
        darray->capacity *= 2;
        darray_resize(darray, darray->capacity);
    }

    darray->elements[darray->length] = element;
    darray->length++;
}

// TODO: Convert to macro
void* darray_get(darray_t* darray, size_t idx) {
    return darray->elements[idx];
}

void darray_clear(darray_t* darray) {
    darray->length = 0;
}

void darray_copy(darray_t* dst, darray_t* src) {
    darray_clear(dst);
    for (size_t i = 0; i < src->length; i++) {
        darray_append(dst, darray_get(src, i));
    }
}

void darray_free(darray_t* darray) {
    free(darray->elements);
    free(darray);
}

#endif

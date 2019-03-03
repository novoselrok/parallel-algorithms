#ifndef _CELL_H
#define _CELL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

#include "common.h"
#include "body.h"

#define N_CELL_CHILDREN 8
#define THETA 0.5
#define AS_CELL_PTR(x) ((cell_t*)(x))

typedef struct _cell_t {
    struct _cell_t* children[N_CELL_CHILDREN];
    body_t* body;
    double mass;
    // Center of mass
    vect3_t cm;
    vect3_t min_bounds;
    vect3_t max_bounds;
    int parent_idx;
    int array_idx;
} cell_t;

typedef struct _mpi_cell_t {
    int parent_idx;
    vect3_t cm;
    vect3_t min_bounds;
    vect3_t max_bounds;
    double mass;
} mpi_cell_t;

extern MPI_Datatype MPI_Cell;

MPI_Datatype MPI_Cell;

void init_MPI_Cell_datatype() {
    int blocklengths[5] = {1, 3, 3, 3, 1};
    MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[5];

    offsets[0] = offsetof(mpi_cell_t, parent_idx);
    offsets[1] = offsetof(mpi_cell_t, cm);
    offsets[2] = offsetof(mpi_cell_t, min_bounds);
    offsets[3] = offsetof(mpi_cell_t, max_bounds);
    offsets[4] = offsetof(mpi_cell_t, mass);

    MPI_Type_create_struct(5, blocklengths, offsets, types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);
}

void init_cell(cell_t* cell) {
    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        cell->children[i] = NULL;
    }
    cell->body = NULL;
    cell->mass = 0.0;
    cell->parent_idx = -1;
    cell->array_idx = -1;
    init_vect(cell->cm);
    init_vect(cell->min_bounds);
    init_vect(cell->max_bounds);
}

mpi_cell_t* pack_cells(darray_t* cells) {
    mpi_cell_t* packed = malloc(sizeof(mpi_cell_t) * cells->length);
    for (size_t i = 0; i < cells->length; i++) {
        cell_t* cell = AS_CELL_PTR(darray_get(cells, i));
        packed[i].parent_idx = cell->parent_idx;
        packed[i].mass = cell->mass;
        memcpy(packed[i].min_bounds, cell->min_bounds, sizeof(vect3_t));
        memcpy(packed[i].max_bounds, cell->max_bounds, sizeof(vect3_t));
        memcpy(packed[i].cm, cell->cm, sizeof(vect3_t));
    }
    return packed;
}

void init_cell_with_bounds(cell_t* cell, vect3_t min_bounds, vect3_t max_bounds) {
    init_cell(cell);
    memcpy(cell->min_bounds, min_bounds, sizeof(vect3_t));
    memcpy(cell->max_bounds, max_bounds, sizeof(vect3_t));
}

void init_cell_bounds_from_parent(cell_t* parent, cell_t* child, int child_idx) {
    // Find the bounds of the cell
    for (int c = 0; c < DIMS; c++) {
        double half_side = (parent->max_bounds[c] - parent->min_bounds[c]) / 2;
        // For each child this generates a separate subblock
        int shift = (child_idx >> c) & 1;
        child->min_bounds[c] = parent->min_bounds[c] + (shift * half_side);
        child->max_bounds[c] = parent->max_bounds[c] - (!shift * half_side);
    }

}

bool cell_contains_position(cell_t* cell, vect3_t pos) {
    if (pos[X] < cell->min_bounds[X] || pos[X] > cell->max_bounds[X]) {
        return false;
    }
    if (pos[Y] < cell->min_bounds[Y] || pos[Y] > cell->max_bounds[Y]) {
        return false;
    }
    if (pos[Z] < cell->min_bounds[Z] || pos[Z] > cell->max_bounds[Z]) {
        return false;
    }
    return true;
}

bool is_external(cell_t* cell) {
    return cell->children[0] == NULL;
}

void construct_empty_tree(cell_t* cell, int current_level, int max_level) {
    if (current_level == max_level) {
        return;
    }

    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        cell_t* child = malloc(sizeof(cell_t));
        cell->children[i] = child;
        init_cell(child);
        init_cell_bounds_from_parent(cell, child, i);
        construct_empty_tree(child, current_level + 1, max_level);
    }
}

void update_empty_cells(cell_t* cell, int current_level, int max_level) {
    if (current_level == max_level) {
        return;
    }

    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        cell_t* child = cell->children[i];
        update_empty_cells(child, current_level + 1, max_level);
    }

    double mass = 0;
    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        cell_t* child = cell->children[i];
        mass += child->mass;

        double combined_mass = cell->mass + child->mass;
        if (combined_mass == 0.0) {
            continue;
        }
        // Update center of mass
        for (int c = 0; c < DIMS; c++) {
            cell->cm[c] = (cell->mass * cell->cm[c] + child->mass * child->cm[c]) / combined_mass;
        }
    }
    cell->mass = mass;
}

void get_leaves(cell_t* cell, cell_t** leaves, int* idx, int current_level, int max_level) {
    if (current_level == max_level - 1) {
        for (int i = 0; i < N_CELL_CHILDREN; i++) {
            leaves[*idx] = cell->children[i];
            (*idx)++;
        }
    } else {
        for (int i = 0; i < N_CELL_CHILDREN; i++) {
            get_leaves(cell->children[i], leaves, idx, current_level + 1, max_level);
        }
    }
}

void free_cell(cell_t* cell) {
    if (cell == NULL) {
        return;
    }

    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        free_cell(cell->children[i]);
    }

    free(cell);
}

bool cell_contains_bounds(cell_t* cell, vect3_t min_bounds, vect3_t max_bounds) {
    for (int c = 0; c < DIMS; c++) {
        if (cell->min_bounds[c] > min_bounds[c] || cell->max_bounds[c] < max_bounds[c]) {
            return false;
        }
    }
    return true;
}

void insert_empty_cell(cell_t *cell, vect3_t min_bounds, vect3_t max_bounds) {
    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        cell_t* child = cell->children[i];
        if (child != NULL && cell_contains_bounds(child, min_bounds, max_bounds)) {
            insert_empty_cell(child, min_bounds, max_bounds);
            return;
        } else if (child == NULL) {
            cell_t* new_child = malloc(sizeof(cell_t));
            cell->children[i] = new_child;
            init_cell_with_bounds(new_child, min_bounds, max_bounds);
            return;
        }
    }
}

void get_cells_to_send(
        cell_t* cell,
        cell_t* parent,
        vect3_t min_bounds,
        vect3_t max_bounds,
        int min_depth,
        int depth,
        darray_t* cells_to_send) {

    if (depth > min_depth) {
        if (depth - min_depth == 1) {
            cell->parent_idx = -1;
        } else {
            cell->parent_idx = parent->array_idx;
        }
        cell->array_idx = (int) cells_to_send->length;
        darray_append(cells_to_send, AS_VOID_PTR(cell));
    }

    vect3_t bounds_center = {(max_bounds[X] - min_bounds[X]) / 2, (max_bounds[Y] - min_bounds[Y]) / 2, (max_bounds[Z] - min_bounds[Z]) / 2};
    double dx = cell->max_bounds[X] - cell->min_bounds[X];
    double dy = cell->max_bounds[Y] - cell->min_bounds[Y];
    double dz = cell->max_bounds[Z] - cell->min_bounds[Z];
    double size = dx + dy + dz;
    double d = distance(cell->cm, bounds_center);

    if (size / d >= THETA) {
        for (int i = 0; i < N_CELL_CHILDREN; i++) {
            if (cell->children[i] != NULL) {
                get_cells_to_send(cell->children[i], cell, min_bounds, max_bounds, min_depth, depth + 1, cells_to_send);
            } else {
                break;
            }
        }
    }
}

darray_t* reconstruct_received_cells(mpi_cell_t* recv_cells, int recv_cells_count) {
    darray_t* all_cells = malloc(sizeof(darray_t));
    darray_t* root_cells = malloc(sizeof(darray_t));
    darray_init(all_cells, (size_t) recv_cells_count);
    darray_init(root_cells, (size_t) recv_cells_count);

    for (int i = 0; i < recv_cells_count; i++) {
        cell_t* cell = malloc(sizeof(cell_t));
        init_cell_with_bounds(cell, recv_cells[i].min_bounds, recv_cells[i].max_bounds);
        cell->mass = recv_cells[i].mass;
        memcpy(cell->cm, recv_cells[i].cm, sizeof(vect3_t));

        if (recv_cells[i].parent_idx == -1) {
            darray_append(root_cells, AS_VOID_PTR(cell));
        } else {
            cell_t* parent = AS_CELL_PTR(darray_get(all_cells, (size_t) recv_cells[i].parent_idx));
            for (int j = 0; j < N_CELL_CHILDREN; j++) {
                if (parent->children[j] == NULL) {
                    parent->children[j] = cell;
                    break;
                }
            }
        }
        darray_append(all_cells, AS_VOID_PTR(cell));
    }
    darray_free(all_cells);
    return root_cells;
}

void insert_cell(cell_t* cell, cell_t* insert) {
    if (cell->mass != 0.0 || insert->mass != 0.0) {
        for (int c = 0; c < DIMS; c++) {
            cell->cm[c] = (cell->mass * cell->cm[c] + insert->mass * insert->cm[c]) / (cell->mass + insert->mass);
        }
        cell->mass += insert->mass;
    }

    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        if (cell->children[i] != NULL && cell_contains_bounds(cell->children[i], insert->min_bounds, insert->max_bounds)) {
            insert_cell(cell->children[i], insert);
            return;
        } else if (cell->children[i] == NULL) {
            cell->children[i] = insert;
            return;
        }
    }
}

void insert_body(cell_t* cell, body_t* body) {
    if (cell->body == NULL && is_external(cell)) {
        // Insert directly
        cell->body = body;
        cell->mass = body->mass;
        memcpy(cell->cm, body->position, sizeof(vect3_t));
    } else {
        // Cell contains a body, split cell
        if (is_external(cell)) {
            // Initialize children
            for (int i = 0; i < N_CELL_CHILDREN; i++) {
                cell_t* child = malloc(sizeof(cell_t));
                cell->children[i] = child;
                init_cell(child);
                init_cell_bounds_from_parent(cell, child, i);
                // Insert old body
                if (cell->body != NULL && cell_contains_position(child, cell->body->position)) {
                    insert_body(child, cell->body);
                    cell->body = NULL;
                }
            }
        }

        // Insert new body
        for (int i = 0; i < N_CELL_CHILDREN; i++) {
            if (cell_contains_position(cell->children[i], body->position)) {
                insert_body(cell->children[i], body);
                break;
            }
        }

        // Update mass and center of mass of the cell
        for (int c = 0; c < DIMS; c++) {
            cell->cm[c] = (cell->mass * cell->cm[c] + body->mass * body->position[c]) / (cell->mass + body->mass);
        }
        cell->mass += body->mass;
    }
}

void add_force(body_t* body, vect3_t position, double mass) {
    double dx = position[X] - body->position[X];
    double dy = position[Y] - body->position[Y];
    double dz = position[Z] - body->position[Z];

    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    double dist_cubed = dist * dist * dist;
    double f = (-G * body->mass * mass) / dist_cubed;

    body->force[X] += f * dx;
    body->force[Y] += f * dy;
    body->force[Z] += f * dz;
}

void compute_force(cell_t* cell, body_t* body) {
    if ((cell->body != NULL && cell->body->id == body->id) || cell->mass == 0.0) {
        return;
    }
    if (is_external(cell) && cell->body != NULL) {
        // Update force acting on body
        add_force(body, cell->body->position, cell->body->mass);
    } else {
        double dx = cell->max_bounds[X] - cell->min_bounds[X];
        double dy = cell->max_bounds[Y] - cell->min_bounds[Y];
        double dz = cell->max_bounds[Z] - cell->min_bounds[Z];
        double size = dx + dy + dz;

        double d = distance(cell->cm, body->position);
        if (size / d < THETA) {
            // Treat this cell as a body
            add_force(body, cell->cm, cell->mass);
        } else {
            // Recursively apply force computation on each child cell
            for (int i = 0; i < N_CELL_CHILDREN; i++) {
                if (cell->children[i] != NULL) {
                    compute_force(cell->children[i], body);
                }
            }
        }
    }
}

#endif

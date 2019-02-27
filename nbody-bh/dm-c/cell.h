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

typedef struct _cell_t {
    struct _cell_t* children[N_CELL_CHILDREN];
    body_t* body;
    double mass;
    // Center of mass
    vect3_t cm;
    vect3_t min_bounds;
    vect3_t max_bounds;
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
    init_vect(cell->cm);
    init_vect(cell->min_bounds);
    init_vect(cell->max_bounds);
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
    if (cell->mass == 0) {
        printf("\n");
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
    if ((cell->body != NULL && cell->body->id == body->id) || cell->mass == 0) {
        return;
    }

    if (is_external(cell)) {
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
                compute_force(cell->children[i], body);
            }
        }
    }
}

#endif

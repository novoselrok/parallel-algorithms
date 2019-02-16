#ifndef _CELL_H
#define _CELL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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

                // Find the bounds of the cell
                for (int c = 0; c < DIMS; c++) {
                    double half_side = (cell->max_bounds[c] - cell->min_bounds[c]) / 2;
                    // For each child this generates a separate subblock
                    int shift = (i >> c) & 1;
                    child->min_bounds[c] = cell->min_bounds[c] + (shift * half_side);
                    child->max_bounds[c] = cell->max_bounds[c] - (!shift * half_side);
                }

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
        for (int c = 0; c < 3; c++) {
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

void free_cell(cell_t* cell) {
    if (cell == NULL) {
        return;
    }

    for (int i = 0; i < N_CELL_CHILDREN; i++) {
        free_cell(cell->children[i]);
    }

    free(cell);
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

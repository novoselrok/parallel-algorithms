#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "common.h"
#include "cell.h"

#define MAX_LINE_LENGTH 1024

void get_universe_size(vect3_t universe_min, vect3_t universe_max, body_t* bodies, int n_bodies) {
    universe_min[0] = universe_min[1] = universe_min[2] = DBL_MAX;
    universe_max[0] = universe_max[1] = universe_max[2] = DBL_MIN;

    for (int i = 0; i < n_bodies; i++) {
        body_t body = bodies[i];
        for (int j = 0; j < DIMS; j++) {
            if (body.position[j] < universe_min[j]) {
                universe_min[j] = body.position[j];
            }
            if (body.position[j] > universe_max[j]) {
                universe_max[j] = body.position[j];
            }
        }
    }
}

int main(int argc, char const *argv[]) {
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int n_bodies = atoi(argv[2]);
    int iterations = atoi(argv[3]);

    // Initialize and read the data structures
    body_t* bodies = (body_t*) malloc(n_bodies * sizeof(body_t));

    char str[MAX_LINE_LENGTH];
    for (int i = 0; i < n_bodies; i++) {
        double x, y, z, vx, vy, vz, mass;
        if (fgets(str, MAX_LINE_LENGTH, f) == NULL) {
            printf("Error reading file.\n");
            exit(1);
        }
        sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz, &mass);
        bodies[i].position[X] = x;
        bodies[i].position[Y] = y;
        bodies[i].position[Z] = z;

        bodies[i].velocity[X] = vx;
        bodies[i].velocity[Y] = vy;
        bodies[i].velocity[Z] = vz;

        bodies[i].mass = mass;
        bodies[i].id = i;
    }

    clock_t begin = clock();
    double dt = 0.1;
    for (int iter = 0; iter < iterations; iter++) {
        vect3_t universe_min, universe_max;
        get_universe_size(universe_min, universe_max, bodies, n_bodies);

        for (int i = 0; i < n_bodies; i++) {
            reset_force(&bodies[i]);
        }

        cell_t* root = malloc(sizeof(cell_t));
        init_cell_with_bounds(root, universe_min, universe_max);

        for (int i = 0; i < n_bodies; i++) {
            insert_body(root, &bodies[i]);
        }

        for (int i = 0; i < n_bodies; i++) {
            compute_force(root, &bodies[i]);
        }

        for (int i = 0; i < n_bodies; i++) {
            body_t* body = &bodies[i];
            body->position[X] += dt * body->velocity[X];
            body->position[Y] += dt * body->velocity[Y];
            body->position[Z] += dt * body->velocity[Z];

            body->velocity[X] += dt / body->mass * body->force[X];
            body->velocity[Y] += dt / body->mass * body->force[Y];
            body->velocity[Z] += dt / body->mass * body->force[Z];
        }

        free_cell(root);
    }
    clock_t end = clock();
    printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    // Output positions and velocities
    FILE* outf;
    outf = fopen("out.txt", "w");
    if(outf == NULL) {
        printf("Could not open output file.\n");
        exit(1);
    }

    for (int i = 0; i < n_bodies; i++) {
        body_t* body = &bodies[i];
        fprintf(outf, "%e %e %e %e %e %e\n", body->position[X], body->position[Y], body->position[Z], body->velocity[X], body->velocity[Y], body->velocity[Z]);
    }
    fclose(outf);
}

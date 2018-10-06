#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X 0
#define Y 1
#define Z 2
#define MAX_LINE_LENGTH 1024
#define G 6.67e-11

// Type definition for an array of three doubles
typedef double vect3_t[3];

int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Not enough arguments: ./main <filename> <iterations> <output file>\n");
        exit(1);
    }

    /* Read initial position, velocity and mass of the bodies */
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int iterations = atoi(argv[2]);

    char str[MAX_LINE_LENGTH];
    fgets(str, MAX_LINE_LENGTH, f);
    int n_bodies;
    sscanf(str, "%d", &n_bodies);

    vect3_t* forces = malloc(n_bodies * sizeof(vect3_t));
    vect3_t* positions = malloc(n_bodies * sizeof(vect3_t));
    vect3_t* velocities = malloc(n_bodies * sizeof(vect3_t));
    double* masses = malloc(n_bodies * sizeof(double));

    for (int i = 0; i < n_bodies; i++) {
        double x, y, z, vx, vy, vz, mass;
        fgets(str, MAX_LINE_LENGTH, f);
        sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz, &mass);
        positions[i][X] = x;
        positions[i][Y] = y;
        positions[i][Z] = z;

        velocities[i][X] = vx;
        velocities[i][Y] = vy;
        velocities[i][Z] = vz;

        masses[i] = mass;
    }

    double dt = 0.1;
    for (int i = 0; i < iterations; i++) {
        // Reset forces
        for (int q = 0; q < n_bodies; q++) {
            forces[q][X] = 0.0;
            forces[q][Y] = 0.0;
            forces[q][Z] = 0.0;
        }
        // Compute forces
        for (int q = 0; q < n_bodies; q++) {
            for (int k = 0; k < n_bodies; k++) {
                // Using Newton's third law of motion we can halve the number of computations needed
                if (k <= q) {
                    continue;
                }
                double x_diff = positions[q][X] - positions[k][X];
                double y_diff = positions[q][Y] - positions[k][Y];
                double z_diff = positions[q][Z] - positions[k][Z];
                double dist = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
                double dist_cubed = dist * dist * dist;
                
                double tmp = -G * masses[q] * masses[k] / dist_cubed;
                double force_qk_x = tmp * x_diff;
                double force_qk_y = tmp * y_diff;
                double force_qk_z = tmp * z_diff;

                forces[q][X] += force_qk_x;
                forces[q][Y] += force_qk_y;
                forces[q][Z] += force_qk_z;
                forces[k][X] -= force_qk_x;
                forces[k][Y] -= force_qk_y;
                forces[k][Z] -= force_qk_z;
            }
        }
        // Update positions and velocities
        for (int q = 0; q < n_bodies; q++) {
            positions[q][X] += dt * velocities[q][X];
            positions[q][Y] += dt * velocities[q][Y];
            positions[q][Z] += dt * velocities[q][Z];
            velocities[q][X] += dt / masses[q] * forces[q][X];
            velocities[q][Y] += dt / masses[q] * forces[q][Y];
            velocities[q][Z] += dt / masses[q] * forces[q][Z];
        }
        // for (int q = 0; q < n_bodies; q++) {
        //     printf("%e %e %e\n", forces[q][X], forces[q][Y], forces[q][Z]);
        // }
        // printf("\n===\n");
    }

    FILE* outf;
    outf = fopen(argv[3], "w");
    if(outf == NULL) {
        printf("Could not open output file.\n");
        exit(1);
    }

    for (int q = 0; q < n_bodies; q++) {
        fprintf(outf, "%e %e %e %e %e %e\n", positions[q][X], positions[q][Y], positions[q][Z], velocities[q][X], velocities[q][Y], velocities[q][Z]);
    }

    // Cleanup
    free(forces);
    free(positions);
    free(velocities);
    free(masses);
    fclose(f);
    fclose(outf);

    return 0;
}

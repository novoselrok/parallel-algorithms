#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define X 0
#define Y 1
#define Z 2
#define MAX_LINE_LENGTH 1024
#define G 6.67e-11

// Type definition for an array of three doubles
typedef double vect3_t[3];

int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Not enough arguments: ./main <filename> <nbodies> <iterations> <output file>\n");
        exit(1);
    }
    int num_threads = 16;
    omp_set_num_threads(num_threads);
    /* Read initial position, velocity and mass of the bodies */
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int n_bodies = atoi(argv[2]);
    int iterations = atoi(argv[3]);

    vect3_t* forces = malloc(n_bodies * sizeof(vect3_t));
    vect3_t* positions = malloc(n_bodies * sizeof(vect3_t));
    vect3_t* velocities = malloc(n_bodies * sizeof(vect3_t));
    double* masses = malloc(n_bodies * sizeof(double));

    char str[MAX_LINE_LENGTH];
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

    vect3_t** thread_loc_forces = malloc(num_threads * sizeof(vect3_t*));
    for (int rank = 0; rank < num_threads; rank++) {
        thread_loc_forces[rank] = malloc(n_bodies * sizeof(vect3_t));
    }
    double dt = 0.1;
    double start_time = omp_get_wtime();
    for (int i = 0; i < iterations; i++) {
        #pragma omp parallel
        {
            int my_rank = omp_get_thread_num();
            // Reset forces
            #pragma omp for
            for (int q = 0; q < n_bodies; q++) {
                forces[q][X] = 0.0;
                forces[q][Y] = 0.0;
                forces[q][Z] = 0.0;
            }
            
            for (int q = 0; q < n_bodies; q++) {
                thread_loc_forces[my_rank][q][X] = 0.0;
                thread_loc_forces[my_rank][q][Y] = 0.0;
                thread_loc_forces[my_rank][q][Z] = 0.0;
            }

            // Compute forces
            #pragma omp for
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

                    thread_loc_forces[my_rank][q][X] += force_qk_x;
                    thread_loc_forces[my_rank][q][Y] += force_qk_y;
                    thread_loc_forces[my_rank][q][Z] += force_qk_z;
                    thread_loc_forces[my_rank][k][X] -= force_qk_x;
                    thread_loc_forces[my_rank][k][Y] -= force_qk_y;
                    thread_loc_forces[my_rank][k][Z] -= force_qk_z;
                }
            }
            // Aggregate all forces
            #pragma omp for
            for (int q = 0; q < n_bodies; q++) {
                for (int rank = 0; rank < num_threads; rank++) {
                    forces[q][X] += thread_loc_forces[rank][q][X];
                    forces[q][Y] += thread_loc_forces[rank][q][Y];
                    forces[q][Z] += thread_loc_forces[rank][q][Z];
                }
            }
            // Update positions and velocities
            #pragma omp for
            for (int q = 0; q < n_bodies; q++) {
                positions[q][X] += dt * velocities[q][X];
                positions[q][Y] += dt * velocities[q][Y];
                positions[q][Z] += dt * velocities[q][Z];
                velocities[q][X] += dt / masses[q] * forces[q][X];
                velocities[q][Y] += dt / masses[q] * forces[q][Y];
                velocities[q][Z] += dt / masses[q] * forces[q][Z];
            }
        }
    }
    printf("time: %.16g\n", omp_get_wtime() - start_time);
    for (int rank = 0; rank < num_threads; rank++) {
        free(thread_loc_forces[rank]);
    }
    free(thread_loc_forces);

    FILE* outf;
    outf = fopen(argv[4], "w");
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

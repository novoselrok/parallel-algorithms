#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char const *argv[]) {
    MPI_Init(NULL, NULL);
    const int n = 100000000;
    double alpha = 0.9;

    int world_size;
    int my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int elements_per_proc = n / world_size;

    int* counts = (int*) malloc(world_size * sizeof(int));
    int* displs = (int*) malloc(world_size * sizeof(int));

    int acc = 0;
    for (int i = 0; i < world_size; i++, acc += elements_per_proc) {
        counts[i] = elements_per_proc;
        displs[i] = acc;
    }

    double* A;
    double* B;
    double* C;
    if (my_rank == 0) {
        A = (double*) malloc(n * sizeof(double));
        B = (double*) malloc(n * sizeof(double));
        C = (double*) malloc(n * sizeof(double));

        srand(42);
        for (int i = 0; i < n; i++) {
            A[i] = ((double) rand()) / RAND_MAX;
            B[i] = ((double) rand()) / RAND_MAX;
        }
    }

    double* my_A = (double*) malloc(elements_per_proc * sizeof(double));
    double* my_B = (double*) malloc(elements_per_proc * sizeof(double));
    double* my_C = (double*) malloc(elements_per_proc * sizeof(double));
    
    double start;
    if (my_rank == 0) {
        start = MPI_Wtime();
    }

    MPI_Scatterv(A, counts, displs, MPI_DOUBLE, my_A, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, counts, displs, MPI_DOUBLE, my_B, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < elements_per_proc; i++) {
        my_C[i] = my_A[i] + alpha * my_B[i];
    }
    MPI_Gatherv(my_C, elements_per_proc, MPI_DOUBLE, C, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Time elapsed: %fs\n", MPI_Wtime() - start);
    }

    MPI_Finalize();
    return 0;
}

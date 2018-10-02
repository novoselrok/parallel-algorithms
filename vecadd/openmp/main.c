#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char const *argv[])
{
    const int n = 100000000;

    double* A = (double*) malloc(n * sizeof(double));
    double* B = (double*) malloc(n * sizeof(double));
    double* C = (double*) malloc(n * sizeof(double));
    double alpha = 0.9;

    srand(42);
    for (int i = 0; i < n; i++) {
        A[i] = ((double) rand()) / RAND_MAX;
        B[i] = ((double) rand()) / RAND_MAX;
    }

    double start = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        C[i] = A[i] + alpha * B[i];
    }
    printf("Time elapsed: %fs\n", omp_get_wtime() - start);

    return 0;
}

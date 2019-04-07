#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1000000
double compute_pi(int n) {
    int withinCircle = 0;
    for(int _i = 0; _i < n; _i++) {
        double x = (rand() / (double) RAND_MAX) * 2 - 1;
        double y = (rand() / (double) RAND_MAX) * 2 - 1;
        double r2 = x * x + y * y;
        if (r2 < 1.0) {
            withinCircle += 1;
        }
    }
    return (withinCircle / (double) n) * 4.0;
}
int main() {
    int nthreads;
    #pragma omp parallel
    {
        #pragma omp master
        nthreads = omp_get_num_threads();
    }
    double result = 0.0;
    #pragma omp parallel for reduction(+: result)
    for(int _i = 0; _i < nthreads; _i++) {
        result = compute_pi((int) ceil(N / nthreads));
    }
    printf("%f\n", result / nthreads);
}


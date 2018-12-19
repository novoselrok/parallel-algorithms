#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int T;

int cmpfunc(const void *a, const void *b) {
   return *(T*)a - *(T*)b;
}

int main(int argc, char const *argv[]) {
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int n = atoi(argv[2]);

    T* arr = malloc(sizeof(T) * n);
    for (int i = 0; i < n; i++) {
        T tmp;
        fscanf(f, "%d", &tmp);
        arr[i] = tmp;
    }
    fclose(f);

    clock_t begin = clock();
    qsort(arr, n, sizeof(T), cmpfunc);
    clock_t end = clock();
	printf("Time elapsed is %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);

    for (int i = 0; i < n - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            printf("Array not sorted!\n");
            exit(1);
        }
    }
    printf("Array sorted!\n");

    return 0;
}

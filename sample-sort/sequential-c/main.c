#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int T;

int cmpfunc(const void *a, const void *b) {
   return *(T*)a - *(T*)b;
}

int partition(T* arr, int left, int right) {
    int pivot = arr[right];
    while (left < right) {
      while (arr[left] < pivot) {
        left++;
      }
      while (arr[right] > pivot) {
        right--;
      }
      if (left <= right) {
        int temp = arr[left];
        arr[left] = arr[right];
        arr[right] = temp;
      }
    }
    return left; // pivot index
}

void myqsort(T* arr, int left, int right) {
    if (left >= right) {
        return;
    }

    int pivot_idx = partition(arr, left, right);
    myqsort(arr, left, pivot_idx - 1);  // sort left of pivot
    myqsort(arr, pivot_idx, right); // sort right of pivot
}

#define REPEAT 100
#define MY_RAND_MAX ((1U << 31) - 1)
int get_random_number(unsigned int seed) {
    return (seed * 1103515245 + 12345) & MY_RAND_MAX;
}

void init_random_array(int* arr, int length, int initial_seed) {
    int random_num = get_random_number(initial_seed);
    for (int i = 0; i < length; i++) {
        arr[i] = random_num;
        random_num = get_random_number(random_num);
    }
}

int main(int argc, char const *argv[]) {
    int n = atoi(argv[1]);
    double times[REPEAT];
    for (int iter = 0; iter < REPEAT; iter++) {
        T* arr = malloc(sizeof(T) * n);
        init_random_array(arr, n, iter + 1); 

        clock_t begin = clock();
        myqsort(arr, 0, n - 1);
        clock_t end = clock();
        times[iter] = (double)(end - begin) / CLOCKS_PER_SEC;

        for (int i = 0; i < n - 1; i++) {
            if (arr[i] > arr[i + 1]) {
                printf("Array not sorted!\n");
                exit(1);
            }
        }
        free(arr);
    }

    double times_sum = 0.0;
    for (int i = 0; i < REPEAT; i++) {
        times_sum += times[i];
    }
    printf("%f\n", times_sum / REPEAT);

    return 0;
}

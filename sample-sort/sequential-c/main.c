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
    // qsort(arr, n, sizeof(T), cmpfunc);
    myqsort(arr, 0, n - 1);
    clock_t end = clock();
	printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);

    for (int i = 0; i < n - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            printf("Array not sorted!\n");
            exit(1);
        }
    }
    return 0;
}

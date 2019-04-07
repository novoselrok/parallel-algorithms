#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

#define OVERSAMPLING_FACTOR 128
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef int T;

int cmpfunc(const void *a, const void *b) {
   return *(T*)a - *(T*)b;
}

/// QuickSort
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
///

int get_number_of_bins() {
    int m;
    #pragma omp parallel 
    {
        #pragma omp master
        m = omp_get_num_threads();    
    }
    return m;
}

void get_sample_keys(T* arr, int n, int m, T* sample_keys) {
    int total_sampled_keys = m * OVERSAMPLING_FACTOR;
    T* sampled_keys = (T*) malloc(sizeof(T) * total_sampled_keys);
    for (int i = 0; i < total_sampled_keys; i++) {
        // int idx = (int) ((rand() / (double) RAND_MAX) * n);
        sampled_keys[i] = arr[i];
    }

    for (int i = 0; i < m - 1; i++) {
        sample_keys[i] = sampled_keys[(i + 1) * OVERSAMPLING_FACTOR];
    }
    // qsort(sample_keys, m - 1, sizeof(T), cmpfunc);
    myqsort(sample_keys, 0, m - 2);

    // Cleanup
    free(sampled_keys);
}

int binary_search(T* arr, int n, T el) {
    int left = 0;
    int right = n;

    while (left < right) {
        int middle = (int) (left + right) / 2;
        if (arr[middle] >= el) {
            right = middle;
        } else {
            left = middle + 1;            
        }
    }

    return left;
}

void map_keys_to_bins(T* arr, int n, int m, T* sample_keys, int* index, int* freq) {
    for (int i = 0; i < n; i++) {
        int bin_idx = binary_search(sample_keys, m - 1, arr[i]);
        index[i] = bin_idx;
        freq[bin_idx] += 1;
    }
}

void bin(T* arr, int n, T*** bins, int** tally, int m) {
    T* sample_keys = (T*) malloc(sizeof(T) * (m - 1));
    get_sample_keys(arr, n, m, sample_keys);

    int block_size = (int) (n + m - 1) / m;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * block_size;
        int end = MIN(start + block_size, n);
        int subarray_size = end - start;
        
        int* freq = (int*) malloc(sizeof(int) * m);
        int* index = (int*) malloc(sizeof(int) * subarray_size);

        for (int i = 0; i < m; i++) {
            freq[i] = 0;
        }

        map_keys_to_bins(&arr[start], subarray_size, m, sample_keys, index, freq);

        // Initialize bins
        for (int i = 0; i < m; i++) {
            bins[thread_id][i] = (T*) malloc(sizeof(T) * freq[i]);
        }

        for (int i = 0; i < subarray_size; i++) {
            int bin_idx = index[i];
            int key = arr[start + i];
            // We use freq[bin_idx] as an index into the bin array
            freq[bin_idx] -= 1;
            bins[thread_id][bin_idx][freq[bin_idx]] = key;
            tally[thread_id][bin_idx] += 1;
        }
        // Cleanup
        free(freq);
        free(index);
    }

    // Cleanup
    free(sample_keys);
}

void subsort(T* sorted_array, T*** bins, int** tally, int m) {
    int* col_sum = (int*) malloc(sizeof(int) * m);

    for (int i = 0; i < m; i++) {
        col_sum[i] = 0;
        for (int j = 0; j < m; j++) {
            col_sum[i] += tally[j][i];
        }
    }

    int* prefix_col_sum = (int*) malloc(sizeof(int) * m);
    prefix_col_sum[0] = 0;
    for (int i = 1; i < m; i++) {
        prefix_col_sum[i] = prefix_col_sum[i - 1] + col_sum[i - 1];
    }

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        T* subarray = (T*) malloc(sizeof(T) * col_sum[thread_id]);

        int offset = 0;
        for (int i = 0; i < m; i++) {
            T* bin = bins[i][thread_id];
            int count = tally[i][thread_id];
            memcpy(&subarray[offset], bin, sizeof(T) * count);
            offset += count;
        }

        myqsort(subarray, 0, col_sum[thread_id] - 1);

        memcpy(&sorted_array[prefix_col_sum[thread_id]], subarray, sizeof(T) * col_sum[thread_id]);
        free(subarray);
    }
    free(col_sum);
    free(prefix_col_sum);
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

        int m = get_number_of_bins();
        double begin = omp_get_wtime();
        // 3D array containing keys (array elements) broken down into bins
        // Each thread owns #thread_id row when binning (1st phase)
        // and #thread_id column when subsorting (2nd phase)
        T*** bins = (T***) malloc(sizeof(T**) * m);
        for (int i = 0; i < m; i++) {
            bins[i] = (T**) malloc(sizeof(T*) * m);
        }

        // 2D array containg the number of elements in each bin
        int** tally = (int**) malloc(sizeof(int*) * m);
        for (int i = 0; i < m; i++) {
            tally[i] = (int*) malloc(sizeof(int) * m);

            for (int j = 0; j < m; j++) {
                tally[i][j] = 0;
            }
        }

        bin(arr, n, bins, tally, m);

        T* sorted_array = (T*) malloc(sizeof(T) * n);
        subsort(sorted_array, bins, tally, m);
        times[iter] = omp_get_wtime() - begin;

        for (int i = 0; i < n - 1; i++) {
            if (sorted_array[i] > sorted_array[i + 1]) {
                printf("Array not sorted!\n");
                exit(1);
            }
        }

        // Cleanup
        free(sorted_array);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                free(bins[i][j]);
            }
            free(tally[i]);
            free(bins[i]);
        }

        free(tally);
        free(bins);
        free(arr);
    }
    double times_sum = 0.0;
    for (int i = 0; i < REPEAT; i++) {
        times_sum += times[i];
    }
    printf("%f\n", times_sum / REPEAT);

    return 0;
}

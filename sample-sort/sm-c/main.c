#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

#define OVERSAMPLING_FACTOR 64
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef int T;

int cmpfunc(const void *a, const void *b) {
   return *(T*)a - *(T*)b;
}

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
        int idx = (int) ((rand() / (double) RAND_MAX) * n);
        sampled_keys[i] = arr[idx];
    }

    for (int i = 0; i < m - 1; i++) {
        sample_keys[i] = sampled_keys[(i + 1) * OVERSAMPLING_FACTOR];
    }
    qsort(sample_keys, m - 1, sizeof(T), cmpfunc);

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
        // printf("Key: %d, Bin idx: %d\n", arr[i], bin_idx);
        index[i] = bin_idx;
        freq[bin_idx] += 1;
    }
}

void bin(T* arr, int n, T*** bins, int** tally, int m) {
    T* sample_keys = (T*) malloc(sizeof(T) * (m - 1));
    get_sample_keys(arr, n, m, sample_keys);
    
    // for (int i = 0; i < m - 1; i++) {
    //     printf("%d ", sample_keys[i]);
    // }
    // printf("\n");
    
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
        // printf("%d ", col_sum[i]);
    }
    // printf("\n");

    int* prefix_col_sum = (int*) malloc(sizeof(int) * m);
    prefix_col_sum[0] = 0;
    for (int i = 1; i < m; i++) {
        prefix_col_sum[i] = prefix_col_sum[i - 1] + col_sum[i - 1];
        // printf("%d ", prefix_col_sum[i]);
    }
    // printf("\n");

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

        qsort(subarray, col_sum[thread_id], sizeof(T), cmpfunc);

        // #pragma omp critical
        // {
        // printf("Thread id: %d\n", thread_id);
        // for (int i = 0; i < col_sum[thread_id]; i++) {
        //     printf("%d ", subarray[i]);
        // }
        // printf("\n");
        // }

        memcpy(&sorted_array[prefix_col_sum[thread_id]], subarray, sizeof(T) * col_sum[thread_id]);
    }
}

int main(int argc, char const *argv[]) {
    srand(time(NULL));

    // T arr[20] = {10, 18, 16, 14, 0, 17, 11, 2, 3, 9, 5, 7, 4, 19, 6, 15, 8, 1, 13, 12};
    // int n = 20;
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

    int m = get_number_of_bins();
    printf("Number of bins: %d\n", m);
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

    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < m; j++) {
    //         T* bin = bins[i][j];
    //         printf("[");
    //         for (int k = 0; k < tally[i][j]; k++) {
    //             printf("%d ", bin[k]);
    //         }
    //         printf("] ");
    //     }
    //     printf("\n");
    // }
    T* sorted_array = (T*) malloc(sizeof(T) * n);
    subsort(sorted_array, bins, tally, m);
    printf("%.16g\n", omp_get_wtime() - begin);

    for (int i = 0; i < n - 1; i++) {
        // printf("%d ", sorted_array[i]);
        if (sorted_array[i] > sorted_array[i + 1]) {
            printf("Array not sorted!\n");
            exit(1);
        }
    }
    // printf("%d\n", sorted_array[n - 1]);
    printf("Array sorted!\n");

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
    return 0;
}

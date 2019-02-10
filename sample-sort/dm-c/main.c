#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <stddef.h>

#define OVERSAMPLING_FACTOR 128
#define MASTER 0
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

void get_sample_keys(T* arr, int n, int m, T* sample_keys) {
    int total_sampled_keys = m * OVERSAMPLING_FACTOR;
    T* sampled_keys = (T*) malloc(sizeof(T) * total_sampled_keys);
    for (int i = 0; i < total_sampled_keys; i++) {
        sampled_keys[i] = arr[i];
    }

    for (int i = 0; i < m - 1; i++) {
        sample_keys[i] = sampled_keys[(i + 1) * OVERSAMPLING_FACTOR];
    }
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

void compute_bins(int* arr, int nkeys, int m, int* sample_keys, int** bins, int* tally) {
    int* index = (int*) malloc(sizeof(int) * nkeys);
    int* freq = (int*) malloc(sizeof(int) * m);
    for (int i = 0; i < m; i++) {
        freq[i] = 0;
    }

    map_keys_to_bins(arr, nkeys, m, sample_keys, index, freq);

    for (int i = 0; i < m; i++) {
        bins[i] = malloc(sizeof(T) * freq[i]);
    }

    for (int i = 0; i < nkeys; i++) {
        int bin_idx = index[i];
        int key = arr[i];
        // We use freq[bin_idx] as an index into the bin array
        freq[bin_idx] -= 1;
        bins[bin_idx][freq[bin_idx]] = key;
        tally[bin_idx] += 1;
    }

    free(freq);
    free(index);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Read data
    srand(time(NULL));
    /* Read the points */
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int nkeys = atoi(argv[2]);
    int m = world_size;

    int* sendcounts = malloc(world_size * sizeof(int));
    int* displs = malloc(world_size * sizeof(int));

    int block_size = (int) (nkeys + world_size - 1) / world_size;
    for (int i = 0; i < world_size; i++) {
        int start = i * block_size;
        int end = MIN(start + block_size, nkeys);
        int subarray_size = end - start;
        
        sendcounts[i] = subarray_size;
        displs[i] = start;
    }

    int* sample_keys = (int*) malloc((m - 1) * sizeof(int));
    int* arr = NULL;
    int* sorted_array = NULL;
    int my_nkeys = sendcounts[rank];
    int* my_arr = (int*) malloc(my_nkeys * sizeof(int));
    if (rank == MASTER) {
        arr = malloc(nkeys * sizeof(int));
        for (int i = 0; i < nkeys; i++) {
            int tmp;
            fscanf(f, "%d", &tmp);
            arr[i] = tmp;
        }
        fclose(f);
        get_sample_keys(arr, nkeys, m, sample_keys);
        sorted_array = (int*) malloc(nkeys * sizeof(int));
    }

    MPI_Scatterv(
        arr, // sendbuf
        sendcounts,
        displs,
        MPI_INT,

        my_arr, // recvbuf
        my_nkeys,
        MPI_INT,
        MASTER,

        MPI_COMM_WORLD
    );

    free(arr);
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Bcast(
        sample_keys,
        m - 1,
        MPI_INT,
        MASTER,
        MPI_COMM_WORLD
    );

    int** my_bins = (int**) malloc(m * sizeof(int*));
    int* my_tally = (int*) malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) {
        my_tally[i] = 0;
    }

    compute_bins(my_arr, my_nkeys, m, sample_keys, my_bins, my_tally);
    int* bin_counts = (int*) malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; i++) {
        if (i == rank) {
            bin_counts[i] = my_tally[i];
        } else {
            MPI_Request req;
            MPI_Isend(
                &my_tally[i],
                1,
                MPI_INT,
                i,
                rank,
                MPI_COMM_WORLD,
                &req
            );
            MPI_Request_free(&req);
        }
    }

    for (int source = 0; source < world_size; source++) {
        if (source == rank) continue;

        int recv_count;
        MPI_Status status;
        MPI_Recv(
            &recv_count,
            1,
            MPI_INT,
            source,
            MPI_ANY_TAG,
            MPI_COMM_WORLD,
            &status
        );
        bin_counts[source] = recv_count;
    }

    int subarray_size = 0;
    for (int i = 0; i < world_size; i++) {
        subarray_size += bin_counts[i];
    }
    int* subarray = (int*) malloc(subarray_size * sizeof(int));
    int offset = 0;

    for (int i = 0; i < world_size; i++) {
        if (i == rank) {
            memcpy(&subarray[offset], my_bins[i], bin_counts[i] * sizeof(int));
            offset += bin_counts[i];
        } else {
            MPI_Request req;
            MPI_Isend(
                my_bins[i],
                my_tally[i],
                MPI_INT,
                i,
                rank,
                MPI_COMM_WORLD,
                &req
            );
            MPI_Request_free(&req);
        }
    }

    for (int source = 0; source < world_size; source++) {
        if (source == rank) continue;

        int* recv_array = (int*) malloc(bin_counts[source] * sizeof(int));
        MPI_Status status;
        MPI_Recv(
            recv_array,
            bin_counts[source],
            MPI_INT,
            source,
            MPI_ANY_TAG,
            MPI_COMM_WORLD,
            &status
        );
        memcpy(&subarray[offset], recv_array, bin_counts[source] * sizeof(int));
        offset += bin_counts[source];
    }

    myqsort(subarray, 0, subarray_size - 1);

    int* subarray_sizes = NULL;
    if (rank != MASTER) {
        MPI_Send(
            &subarray_size,
            1,
            MPI_INT,
            MASTER,
            rank,
            MPI_COMM_WORLD
        );
    } else {
        subarray_sizes = (int*) malloc(world_size * sizeof(int));
        subarray_sizes[MASTER] = subarray_size;
        for (int source = 0; source < world_size; source++) {
            if (source == MASTER) continue;
            int source_subarray_size;
            MPI_Status status;
            MPI_Recv(
                &source_subarray_size,
                1,
                MPI_INT,
                source,
                MPI_ANY_TAG,
                MPI_COMM_WORLD,
                &status
            );
            subarray_sizes[source] = source_subarray_size;
        }
    }

    if (rank != MASTER) {
        MPI_Send(
            subarray,
            subarray_size,
            MPI_INT,
            MASTER,
            rank,
            MPI_COMM_WORLD
        );
    } else {
        int offset = 0;
        memcpy(&sorted_array[offset], subarray, subarray_size * sizeof(int));
        offset += subarray_size;

        for (int source = 0; source < world_size; source++) {
            if (source == MASTER) continue;
            int* source_subarray = (int*) malloc(subarray_sizes[source] * sizeof(int));
            MPI_Status status;
            MPI_Recv(
                source_subarray,
                subarray_sizes[source],
                MPI_INT,
                source,
                MPI_ANY_TAG,
                MPI_COMM_WORLD,
                &status
            );
            memcpy(&sorted_array[offset], source_subarray, subarray_sizes[source] * sizeof(int));
            offset += subarray_sizes[source];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == MASTER) {
        printf("%f\n", end_time - start_time);
        for (int i = 0; i < nkeys - 1; i++) {
            if (sorted_array[i] > sorted_array[i + 1]) {
                printf("Array not sorted!\n");
                exit(1);
            }
        }
    }

    MPI_Finalize();
}

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <stddef.h>

#define POINT_SIZE 100
#define MASTER 0
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

// Type definition for a point of size POINT_SIZE
typedef double point_t[POINT_SIZE];

typedef struct _cluster_t {
    point_t sum;
    point_t mean;
    int count;
} cluster_t;

void create_cluster_mpi_datatype(MPI_Datatype* cluster_datatype) {
    cluster_t tmp_cluster;
    MPI_Datatype tmp_type;

    int blocklen[3] = { POINT_SIZE, POINT_SIZE, 1 };
    MPI_Datatype type[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
    MPI_Aint disp[3];
    disp[0] = offsetof(cluster_t, sum);
    disp[1] = offsetof(cluster_t, mean);
    disp[2] = offsetof(cluster_t, count);
    MPI_Type_create_struct(3, blocklen, disp, type, &tmp_type);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(tmp_type, &lb, &extent);
    MPI_Type_create_resized(tmp_type, lb, extent, cluster_datatype);
    MPI_Type_commit(cluster_datatype);
}

void set_cluster_mean(cluster_t* cluster, point_t mean) {
    memcpy(cluster->mean, mean, POINT_SIZE * sizeof(double));
}

void add_point(cluster_t* cluster, point_t point) {
    for (int i = 0; i < POINT_SIZE; i++) {
        cluster->sum[i] += point[i];
    }
    cluster->count++;
}

double distance(cluster_t* cluster, point_t point) {
    double sum_squares = 0.0;
    for (int i = 0; i < POINT_SIZE; i++) {
        double diff = cluster->mean[i] - point[i];
        sum_squares += (diff * diff);
    }
    return sqrt(sum_squares);
}

void calculate_means(cluster_t* clusters, cluster_t* new_clusters, int k) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < POINT_SIZE; j++) {
            clusters[i].mean[j] = new_clusters[i].sum[j] / new_clusters[i].count;
        }
    }
}

void reset_clusters(cluster_t* clusters, int k) {
    for (int i = 0; i < k; i++) {
        clusters[i].count = 0;
        for (int j = 0; j < POINT_SIZE; j++) {
            clusters[i].sum[j] = 0.0;
        }
    }
}

void print_labels(char* filename, int* labels, int length) {
    FILE* outf;
    outf = fopen(filename, "w");
    if(outf == NULL) {
        printf("Could not open output file.\n");
        exit(1);
    }
    for (int i = 0; i < length; i++) {
        fprintf(outf, "%d ", labels[i]);
    }
    fprintf(outf, "\n");
    fclose(outf);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype cluster_datatype;
    create_cluster_mpi_datatype(&cluster_datatype);

    // Read data
    srand(time(NULL));
    /* Read the points */
    FILE* f;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        printf("File does not exist.\n");
        exit(1);
    }
    int k = atoi(argv[2]);
    int max_iter = atoi(argv[3]);
    int n_points = atoi(argv[4]);

    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    int* sendcounts = malloc(world_size * sizeof(int));
    int* displs = malloc(world_size * sizeof(int));
    int* recvcounts = malloc(world_size * sizeof(int));
    int* recvdispls = malloc(world_size * sizeof(int));
    
    int block_size = (int) (n_points + world_size - 1) / world_size;
    for (int i = 0; i < world_size; i++) {
        int start = i * block_size;
        int end = MIN(start + block_size, n_points);
        int subarray_size = end - start;
        
        sendcounts[i] = subarray_size * POINT_SIZE;
        displs[i] = start * POINT_SIZE;
        recvcounts[i] = subarray_size;
        recvdispls[i] = start;
    }

    point_t* points = NULL;
    cluster_t* clusters = malloc(k * sizeof(cluster_t));
    int my_points_count = recvcounts[rank];
    point_t* my_points = malloc(my_points_count * sizeof(point_t));
    int* my_labels = malloc(my_points_count * sizeof(int));

    if (rank == MASTER) {
        points = malloc(n_points * sizeof(point_t));
        for (int i = 0; i < n_points; i++) {
            for (int j = 0; j < POINT_SIZE; j++) {
                double tmp;
                fscanf(f, "%lf", &tmp);
                points[i][j] = tmp;
            }
        }
        fclose(f);

        for (int i = 0; i < k; i++) {
            int idx = (int) ((rand() / (double) RAND_MAX) * n_points);
            set_cluster_mean(&clusters[i], points[idx]);
        }
    }

    MPI_Scatterv(
        points, // sendbuf
        sendcounts,
        displs,
        MPI_DOUBLE,

        my_points, // recvbuf
        sendcounts[rank],
        MPI_DOUBLE,
        MASTER,

        MPI_COMM_WORLD
    );

    MPI_Bcast(
        clusters,
        k,
        cluster_datatype,
        MASTER,
        MPI_COMM_WORLD
    );

    cluster_t* new_clusters = malloc(k * sizeof(cluster_t));
    cluster_t* gathered_clusters = NULL;
    if (rank == MASTER) {
        gathered_clusters = malloc(k * world_size * sizeof(cluster_t));
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // (Re)set new_clusters before use
        reset_clusters(new_clusters, k);

        for (int i = 0; i < my_points_count; i++) {
            double min_distance = DBL_MAX;
            int min_index = 0;
            // Calculate closest cluster
            for (int j = 0; j < k; j++) {
                double distance_to_cluster = distance(&clusters[j], my_points[i]);
                if (distance_to_cluster < min_distance) {
                    min_distance = distance_to_cluster;
                    min_index = j;
                }
            }
            // Get min distance cluster
            my_labels[i] = min_index;
            add_point(&new_clusters[min_index], my_points[i]);
        }

        MPI_Gather(
            new_clusters,
            k,
            cluster_datatype,

            gathered_clusters,
            k,
            cluster_datatype,

            MASTER,
            MPI_COMM_WORLD
        );

        if (rank == MASTER) {
            for (int i = 1; i < world_size; i++) {
                for (int j = 0; j < k; j++) {
                    gathered_clusters[j].count += gathered_clusters[i * k + j].count;

                    for (int l = 0; l < POINT_SIZE; l++) {
                        gathered_clusters[j].sum[l] += gathered_clusters[i * k + j].sum[l];
                    }
                }
            }

            calculate_means(clusters, gathered_clusters, k);
        }

        MPI_Bcast(
            clusters,
            k,
            cluster_datatype,
            MASTER,
            MPI_COMM_WORLD
        );
    }

    int* labels = NULL;
    if (rank == MASTER) {
        labels = malloc(n_points * sizeof(int));
    }

    MPI_Gatherv(
        my_labels,
        my_points_count,
        MPI_INT,
        
        labels,
        recvcounts,
        recvdispls,
        MPI_INT,

        MASTER,
        MPI_COMM_WORLD
    );

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == MASTER) {
        printf("%f\n", end_time - start_time);
        print_labels("out.txt", labels, n_points);
    }

    if (rank == MASTER) {
        free(gathered_clusters);
        free(points);
        free(labels);
    }

    free(my_labels);
    free(clusters);
    free(new_clusters);
    free(my_points);
    free(sendcounts);
    free(displs);
    free(recvcounts);
    free(recvdispls);

   MPI_Finalize();
}

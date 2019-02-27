#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

#define VERBOSE 1
#define POINT_SIZE 10

// Type definition for a point of size POINT_SIZE
typedef double point_t[POINT_SIZE];

typedef struct _cluster_t {
    point_t sum;
    point_t mean;
    int count;
} cluster_t;

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

void reset_cluster(cluster_t* cluster) {
    cluster->count = 0;
    for (int i = 0; i < POINT_SIZE; i++) {
        cluster->sum[i] = 0.0;
    }
}

void reset_clusters(cluster_t* clusters, int k) {
    for (int i = 0; i < k; i++) {
        reset_cluster(&clusters[i]);
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

int main(int argc, char const *argv[]) {
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

    point_t* points = malloc(n_points * sizeof(point_t));
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < POINT_SIZE; j++) {
            double tmp;
            fscanf(f, "%lf", &tmp);
            points[i][j] = tmp;
        }
    }
    fclose(f);

    /* Algorithm */
    int* labels = malloc(n_points * sizeof(int));

    cluster_t* clusters = malloc(k * sizeof(cluster_t));
    cluster_t* new_clusters = malloc(k * sizeof(cluster_t));
    reset_clusters(clusters, k);
    reset_clusters(new_clusters, k);

    clock_t begin = clock();

    for (int i = 0; i < k; i++) {
        int idx = (int) ((rand() / (double) RAND_MAX) * n_points);
        set_cluster_mean(&clusters[i], points[idx]);
    }

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n_points; i++) {
            double min_distance = DBL_MAX;
            int min_index = 0;
            // Calculate closest cluster
            for (int j = 0; j < k; j++) {
                double distance_to_cluster = distance(&clusters[j], points[i]);
                if (distance_to_cluster < min_distance) {
                    min_distance = distance_to_cluster;
                    min_index = j;
                }
            }
            // Get min distance cluster
            labels[i] = min_index;
            add_point(&new_clusters[min_index], points[i]);
        }

        // Calculate means
        calculate_means(clusters, new_clusters, k);

        // Reset before reusing
        reset_clusters(new_clusters, k);
    }

    clock_t end = clock();
	printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    
    if (VERBOSE) {
        print_labels("out.txt", labels, n_points);
    }

    // Cleanup
    free(points);
    free(labels);
    free(clusters);
    free(new_clusters);
}

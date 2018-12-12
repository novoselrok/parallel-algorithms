#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define VERBOSE 1
#define POINT_SIZE 100

// Type definition for a point of size POINT_SIZE
typedef double point_t[POINT_SIZE];

/**
 * Sums up two points and modifies p1 in place.
 */
void add_point(point_t p1, point_t p2) {
    for (int i = 0; i < POINT_SIZE; i++) {
        p1[i] += p2[i];
    }
}

/**
 * Calculates euclidian distance of two points
 */
double distance(point_t p1, point_t p2) {
    double sum_squares = 0.0;
    for (int i = 0; i < POINT_SIZE; i++) {
        double diff = p1[i] - p2[i];
        sum_squares += (diff * diff);
    }
    return sqrt(sum_squares);
}

/**
 * Returns the index of the minimum value in the array.
 */
int argmin(double* arr, int length) {
    int min_index = 0;
    int min_value = arr[0];
    for (int i = 1; i < length; i++) {
        if (arr[i] < min_value) {
            min_value = arr[i];
            min_index = i;
        }
    }
    return min_index;
}

/**
 * Resets the provided sums and counts.
 */
void reset_cluster_centers_sum(point_t* cluster_centers_sum, int* cluster_centers_count, int k) {
    for (int i = 0; i < k; i++) {
        cluster_centers_count[i] = 0;
        for (int j = 0; j < POINT_SIZE; j++) {
            cluster_centers_sum[i][j] = 0.0;
        }
    }
}

/**
 * Calculates the mean from sums and counts.
 */
void calculate_means(point_t* cluster_centers, point_t* cluster_centers_sum, int* cluster_centers_count, int k) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < POINT_SIZE; j++) {
            cluster_centers[i][j] = cluster_centers_sum[i][j] / cluster_centers_count[i];
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
    point_t* cluster_centers = malloc(k * sizeof(point_t));
    point_t* cluster_centers_sum = malloc(k * sizeof(point_t));
    int* cluster_centers_count = malloc(k * sizeof(int));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < POINT_SIZE; j++) {
            cluster_centers[i][j] = 0.0;
        }
    }

    clock_t begin = clock();
    // Randomly initialize the cluster centers
    for (int i = 0; i < k; i++) {
        int idx = (int) ((rand() / (double) RAND_MAX) * n_points);
        add_point(cluster_centers[i], points[idx]);
    }

    double* cluster_distances = malloc(k * sizeof(double));
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n_points; i++) {
            // Calculate closest cluster
            for (int j = 0; j < k; j++) {
                cluster_distances[j] = distance(points[i], cluster_centers[j]);
            }
            // Get min distance cluster
            labels[i] = argmin(cluster_distances, k);
            add_point(cluster_centers_sum[labels[i]], points[i]);
            cluster_centers_count[labels[i]]++;
        }

        // Calculate means
        calculate_means(cluster_centers, cluster_centers_sum, cluster_centers_count, k);

        // Reset before reusing
        reset_cluster_centers_sum(cluster_centers_sum, cluster_centers_count, k);
    }

    clock_t end = clock();
	printf("Time elapsed is %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
    
    if (VERBOSE) {
        print_labels("out.txt", labels, n_points);
    }

    // Cleanup
    free(points);
    free(labels);
    free(cluster_centers_count);
    free(cluster_centers);
    free(cluster_distances);
}

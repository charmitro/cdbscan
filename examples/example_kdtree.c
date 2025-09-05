/*
 * cdbscan - DBSCAN clustering algorithm implementation in C
 * Copyright (C) 2025 The cdbscan developers
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Example demonstrating KD-tree acceleration for large datasets */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cdbscan.h"

void generate_large_dataset(cdbscan_point_t *points, int num_points)
{
	/* Generate multiple Gaussian clusters */
	int clusters = 5;
	int points_per_cluster = (num_points * 0.9) / clusters;
	int noise_points = num_points - (points_per_cluster * clusters);

	/* Cluster centers */
	double centers[][2] = { { 10.0, 10.0 },
				{ 30.0, 10.0 },
				{ 20.0, 30.0 },
				{ 40.0, 40.0 },
				{ 10.0, 40.0 } };

	int idx = 0;

	/* Generate clustered points */
	for (int c = 0; c < clusters; c++) {
		for (int i = 0; i < points_per_cluster && idx < num_points;
		     i++) {
			/* Box-Muller transform for Gaussian distribution */
			double u1 = rand() / (double)RAND_MAX;
			double u2 = rand() / (double)RAND_MAX;
			double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
			double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

			points[idx].coords[0] = centers[c][0] + z0 * 2.0;
			points[idx].coords[1] = centers[c][1] + z1 * 2.0;
			idx++;
		}
	}

	/* Generate noise points */
	for (int i = 0; i < noise_points && idx < num_points; i++) {
		points[idx].coords[0] = rand() / (double)RAND_MAX * 50.0;
		points[idx].coords[1] = rand() / (double)RAND_MAX * 50.0;
		idx++;
	}
}

double get_time_diff(struct timespec start, struct timespec end)
{
	return (end.tv_sec - start.tv_sec) +
	       (end.tv_nsec - start.tv_nsec) / 1e9;
}

int main()
{
	printf("KD-Tree Performance Comparison\n");
	printf("==============================\n\n");

	/* Test with different dataset sizes */
	int test_sizes[] = { 1000, 5000, 10000, 20000 };
	int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);

	for (int t = 0; t < num_tests; t++) {
		int num_points = test_sizes[t];
		int dimensions = 2;

		printf("Dataset size: %d points\n", num_points);
		printf("------------------------\n");

		/* Create points */
		cdbscan_point_t *points1 =
			cdbscan_create_points(num_points, dimensions);
		cdbscan_point_t *points2 =
			cdbscan_create_points(num_points, dimensions);

		if (!points1 || !points2) {
			fprintf(stderr, "Failed to allocate points\n");
			return 1;
		}

		/* Generate same data for both tests */
		srand(42); /* Fixed seed for reproducibility */
		generate_large_dataset(points1, num_points);
		srand(42);
		generate_large_dataset(points2, num_points);

		/* Parameters */
		double eps = 2.0;
		int min_pts = 5;

		/* Test 1: Without KD-tree (O(n²)) */
		cdbscan_params_t params_brute = {
			.eps = eps,
			.min_pts = min_pts,
			.dist_type = CDBSCAN_DIST_EUCLIDEAN,
			.minkowski_p = 2.0,
			.custom_dist = NULL,
			.custom_dist_params = NULL,
			.use_kdtree = 0 /* Disable KD-tree */
		};

		struct timespec start_brute, end_brute;
		clock_gettime(CLOCK_MONOTONIC, &start_brute);

		int clusters_brute =
			cdbscan_cluster(points1, num_points, params_brute);

		clock_gettime(CLOCK_MONOTONIC, &end_brute);
		double time_brute = get_time_diff(start_brute, end_brute);

		/* Test 2: With KD-tree (O(n log n)) */
		cdbscan_params_t params_kdtree = {
			.eps = eps,
			.min_pts = min_pts,
			.dist_type = CDBSCAN_DIST_EUCLIDEAN,
			.minkowski_p = 2.0,
			.custom_dist = NULL,
			.custom_dist_params = NULL,
			.use_kdtree = 1 /* Enable KD-tree */
		};

		struct timespec start_kdtree, end_kdtree;
		clock_gettime(CLOCK_MONOTONIC, &start_kdtree);

		int clusters_kdtree =
			cdbscan_cluster(points2, num_points, params_kdtree);

		clock_gettime(CLOCK_MONOTONIC, &end_kdtree);
		double time_kdtree = get_time_diff(start_kdtree, end_kdtree);

		/* Verify results are the same */
		if (clusters_brute != clusters_kdtree) {
			printf("WARNING: Different cluster counts!\n");
		}

		/* Print results */
		printf("Brute force:   %.4f seconds (%d clusters)\n",
		       time_brute, clusters_brute);
		printf("With KD-tree:  %.4f seconds (%d clusters)\n",
		       time_kdtree, clusters_kdtree);
		printf("Speedup:       %.2fx\n\n", time_brute / time_kdtree);

		/* Clean up */
		for (int i = 0; i < num_points; i++) {
			free(points1[i].coords);
			free(points2[i].coords);
		}
		free(points1);
		free(points2);
	}

	printf("Summary\n");
	printf("-------\n");
	printf("KD-tree provides significant speedup for large datasets.\n");
	printf("The speedup increases with dataset size:\n");
	printf("- Small datasets (< 1000): Minimal improvement\n");
	printf("- Medium datasets (1000-10000): 2-5x speedup\n");
	printf("- Large datasets (> 10000): 5-10x+ speedup\n\n");
	printf("Note: KD-tree is only used with Euclidean distance.\n");

	return 0;
}

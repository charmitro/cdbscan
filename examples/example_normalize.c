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

/* Example demonstrating data normalization */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdbscan.h"

void generate_unscaled_data(cdbscan_point_t *points, int num_points)
{
	/* Generate clusters with very different scales */
	for (int i = 0; i < num_points; i++) {
		if (i < num_points / 3) {
			/* Cluster 1: small scale (0-1 range) */
			points[i].coords[0] =
				0.5 + (rand() / (double)RAND_MAX - 0.5) * 0.3;
			points[i].coords[1] =
				0.5 + (rand() / (double)RAND_MAX - 0.5) * 0.3;
		} else if (i < 2 * num_points / 3) {
			/* Cluster 2: large scale (100-200 range) */
			points[i].coords[0] =
				150.0 +
				(rand() / (double)RAND_MAX - 0.5) * 30.0;
			points[i].coords[1] =
				150.0 +
				(rand() / (double)RAND_MAX - 0.5) * 30.0;
		} else {
			/* Cluster 3: medium scale (10-20 range) */
			points[i].coords[0] =
				15.0 + (rand() / (double)RAND_MAX - 0.5) * 3.0;
			points[i].coords[1] =
				15.0 + (rand() / (double)RAND_MAX - 0.5) * 3.0;
		}
	}
}

void print_data_stats(cdbscan_point_t *points, int num_points,
		      const char *label)
{
	double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
	double sum_x = 0, sum_y = 0;

	for (int i = 0; i < num_points; i++) {
		if (points[i].coords[0] < min_x)
			min_x = points[i].coords[0];
		if (points[i].coords[0] > max_x)
			max_x = points[i].coords[0];
		if (points[i].coords[1] < min_y)
			min_y = points[i].coords[1];
		if (points[i].coords[1] > max_y)
			max_y = points[i].coords[1];
		sum_x += points[i].coords[0];
		sum_y += points[i].coords[1];
	}

	printf("%s:\n", label);
	printf("  X range: [%.2f, %.2f], mean: %.2f\n", min_x, max_x,
	       sum_x / num_points);
	printf("  Y range: [%.2f, %.2f], mean: %.2f\n", min_y, max_y,
	       sum_y / num_points);
}

int main()
{
	int num_points = 90;
	int dimensions = 2;

	printf("Data Normalization Example\n");
	printf("==========================\n\n");

	/* Test 1: Without normalization */
	printf("Test 1: Without Normalization\n");
	printf("-----------------------------\n");

	cdbscan_point_t *points1 =
		cdbscan_create_points(num_points, dimensions);
	if (!points1) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	generate_unscaled_data(points1, num_points);
	print_data_stats(points1, num_points, "Original data");

	cdbscan_params_t params = {
		.eps = 30.0, /* Large eps needed for unscaled data */
		.min_pts = 4,
		.dist_type = CDBSCAN_DIST_EUCLIDEAN,
		.minkowski_p = 2.0,
		.custom_dist = NULL,
		.custom_dist_params = NULL
	};

	int clusters1 = cdbscan_cluster(points1, num_points, params);
	printf("Clusters found: %d (eps=%.1f)\n\n", clusters1, params.eps);

	/* Test 2: With min-max normalization */
	printf("Test 2: With Min-Max Normalization\n");
	printf("-----------------------------------\n");

	cdbscan_point_t *points2 =
		cdbscan_create_points(num_points, dimensions);
	if (!points2) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	generate_unscaled_data(points2, num_points);
	cdbscan_normalize_minmax(points2, num_points);
	print_data_stats(points2, num_points, "After min-max normalization");

	params.eps = 0.3; /* Much smaller eps works after normalization */
	int clusters2 = cdbscan_cluster(points2, num_points, params);
	printf("Clusters found: %d (eps=%.1f)\n\n", clusters2, params.eps);

	/* Test 3: With z-score normalization */
	printf("Test 3: With Z-Score Normalization\n");
	printf("-----------------------------------\n");

	cdbscan_point_t *points3 =
		cdbscan_create_points(num_points, dimensions);
	if (!points3) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	generate_unscaled_data(points3, num_points);
	cdbscan_normalize_zscore(points3, num_points);
	print_data_stats(points3, num_points, "After z-score normalization");

	params.eps = 1.0; /* Different eps for z-score normalized data */
	int clusters3 = cdbscan_cluster(points3, num_points, params);
	printf("Clusters found: %d (eps=%.1f)\n", clusters3, params.eps);

	/* Clean up */
	for (int i = 0; i < num_points; i++) {
		free(points1[i].coords);
		free(points2[i].coords);
		free(points3[i].coords);
	}
	free(points1);
	free(points2);
	free(points3);

	return 0;
}

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

/* Example demonstrating automatic eps parameter estimation */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdbscan.h"

void generate_test_data(cdbscan_point_t *points, int num_points)
{
	/* Generate well-separated clusters for testing */
	for (int i = 0; i < num_points; i++) {
		if (i < num_points / 3) {
			/* Dense cluster 1 */
			double angle = (double)i / (num_points / 3) * 2 * M_PI;
			points[i].coords[0] =
				2.0 + 0.3 * cos(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.1;
			points[i].coords[1] =
				2.0 + 0.3 * sin(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.1;
		} else if (i < 2 * num_points / 3) {
			/* Dense cluster 2 */
			double angle = (double)(i - num_points / 3) /
				       (num_points / 3) * 2 * M_PI;
			points[i].coords[0] =
				5.0 + 0.3 * cos(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.1;
			points[i].coords[1] =
				2.0 + 0.3 * sin(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.1;
		} else {
			/* Sparse noise */
			points[i].coords[0] = rand() / (double)RAND_MAX * 7.0;
			points[i].coords[1] = rand() / (double)RAND_MAX * 4.0;
		}
	}
}

int main()
{
	int num_points = 150;
	int dimensions = 2;

	printf("Automatic Eps Estimation Example\n");
	printf("================================\n\n");

	/* Create and populate points */
	cdbscan_point_t *points = cdbscan_create_points(num_points, dimensions);
	if (!points) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	generate_test_data(points, num_points);

	/* Estimate eps using k-dist method */
	printf("Step 1: Estimating eps parameter\n");
	printf("---------------------------------\n");

	int k = 4; /* MinPts value - typical for 2D data */
	cdbscan_kdist_result_t *kdist =
		cdbscan_estimate_eps(points, num_points, k);

	if (!kdist) {
		fprintf(stderr, "Failed to estimate eps\n");
		return 1;
	}

	printf("K-value used: %d\n", kdist->k);
	printf("Suggested eps: %.3f\n\n", kdist->suggested_eps);

	/* Show some k-distances for insight */
	printf("Sample k-distances (sorted):\n");
	for (int i = 0; i < 10 && i < num_points; i++) {
		printf("  Point %d: k-dist = %.3f\n", i, kdist->distances[i]);
	}
	printf("  ...\n");
	for (int i = num_points - 5; i < num_points; i++) {
		printf("  Point %d: k-dist = %.3f\n", i, kdist->distances[i]);
	}
	printf("\n");

	/* Test clustering with manual eps */
	printf("Step 2: Testing with manual eps\n");
	printf("--------------------------------\n");

	cdbscan_params_t params_manual = { .eps = 0.2, /* Too small */
					   .min_pts = k,
					   .dist_type = CDBSCAN_DIST_EUCLIDEAN,
					   .minkowski_p = 2.0,
					   .custom_dist = NULL,
					   .custom_dist_params = NULL };

	int clusters_manual =
		cdbscan_cluster(points, num_points, params_manual);

	int noise_manual = 0;
	for (int i = 0; i < num_points; i++) {
		if (points[i].cluster_id == CDBSCAN_NOISE)
			noise_manual++;
	}

	printf("Manual eps = %.3f:\n", params_manual.eps);
	printf("  Clusters: %d, Noise: %d\n\n", clusters_manual, noise_manual);

	/* Reset cluster assignments */
	for (int i = 0; i < num_points; i++) {
		points[i].cluster_id = CDBSCAN_UNCLASSIFIED;
	}

	/* Test clustering with estimated eps */
	printf("Step 3: Testing with estimated eps\n");
	printf("-----------------------------------\n");

	cdbscan_params_t params_auto = { .eps = kdist->suggested_eps,
					 .min_pts = k,
					 .dist_type = CDBSCAN_DIST_EUCLIDEAN,
					 .minkowski_p = 2.0,
					 .custom_dist = NULL,
					 .custom_dist_params = NULL };

	int clusters_auto = cdbscan_cluster(points, num_points, params_auto);

	int noise_auto = 0;
	for (int i = 0; i < num_points; i++) {
		if (points[i].cluster_id == CDBSCAN_NOISE)
			noise_auto++;
	}

	printf("Estimated eps = %.3f:\n", params_auto.eps);
	printf("  Clusters: %d, Noise: %d\n", clusters_auto, noise_auto);
	printf("\nThe estimated eps typically finds the major clusters\n");
	printf("while treating outliers as noise.\n");

	/* Clean up */
	cdbscan_free_kdist_result(kdist);
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);

	return 0;
}

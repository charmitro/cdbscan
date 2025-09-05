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

/* Example demonstrating different distance metrics */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdbscan.h"

void generate_rectangular_clusters(cdbscan_point_t *points, int num_points)
{
	/* Generate data better suited for Manhattan distance */
	for (int i = 0; i < num_points; i++) {
		if (i < num_points / 3) {
			/* Grid-like cluster 1 */
			points[i].coords[0] = 1.0 + (i % 5) * 0.3;
			points[i].coords[1] = 1.0 + (i / 5) * 0.3;
		} else if (i < 2 * num_points / 3) {
			/* Grid-like cluster 2 */
			int idx = i - num_points / 3;
			points[i].coords[0] = 5.0 + (idx % 5) * 0.3;
			points[i].coords[1] = 1.0 + (idx / 5) * 0.3;
		} else {
			/* Random noise */
			points[i].coords[0] = rand() / (double)RAND_MAX * 8.0;
			points[i].coords[1] = rand() / (double)RAND_MAX * 4.0;
		}
	}
}

int main()
{
	int num_points = 150;
	int dimensions = 2;

	/* Create points */
	cdbscan_point_t *points = cdbscan_create_points(num_points, dimensions);
	if (!points) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	/* Generate data */
	generate_rectangular_clusters(points, num_points);

	/* Test with different distance metrics */
	const char *distance_names[] = { "Euclidean", "Manhattan",
					 "Minkowski(p=3)" };
	cdbscan_dist_type_t distance_types[] = { CDBSCAN_DIST_EUCLIDEAN,
						 CDBSCAN_DIST_MANHATTAN,
						 CDBSCAN_DIST_MINKOWSKI };

	printf("Testing Different Distance Metrics\n");
	printf("==================================\n\n");

	for (int d = 0; d < 3; d++) {
		/* Reset cluster assignments */
		for (int i = 0; i < num_points; i++) {
			points[i].cluster_id = CDBSCAN_UNCLASSIFIED;
		}

		/* Set up parameters */
		cdbscan_params_t params = {
			.eps = (d == 1) ? 1.0 :
					  0.8, /* Manhattan needs larger eps */
			.min_pts = 4,
			.dist_type = distance_types[d],
			.minkowski_p = 3.0, /* Used only for Minkowski */
			.custom_dist = NULL,
			.custom_dist_params = NULL
		};

		printf("Distance Metric: %s\n", distance_names[d]);
		printf("Eps: %.2f, MinPts: %d\n", params.eps, params.min_pts);

		/* Run clustering */
		int num_clusters = cdbscan_cluster(points, num_points, params);

		if (num_clusters < 0) {
			fprintf(stderr, "Clustering failed\n");
			continue;
		}

		/* Count results */
		int noise_count = 0;
		for (int i = 0; i < num_points; i++) {
			if (points[i].cluster_id == CDBSCAN_NOISE) {
				noise_count++;
			}
		}

		printf("Clusters found: %d\n", num_clusters);
		printf("Noise points: %d\n", noise_count);
		printf("------------------------\n\n");
	}

	/* Clean up */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);

	return 0;
}

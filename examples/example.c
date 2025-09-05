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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdbscan.h"

/* Generate sample 2D data with clusters */
void generate_sample_data(cdbscan_point_t *points, int num_points)
{
	int cluster_size = num_points / 4;

	/* Generate 3 circular clusters and some noise */
	for (int i = 0; i < num_points; i++) {
		if (i < cluster_size) {
			/* Cluster 1: centered at (2, 2) */
			double angle = (double)i / cluster_size * 2 * M_PI;
			points[i].coords[0] =
				2.0 + 0.5 * cos(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
			points[i].coords[1] =
				2.0 + 0.5 * sin(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
		} else if (i < 2 * cluster_size) {
			/* Cluster 2: centered at (5, 2) */
			double angle = (double)(i - cluster_size) /
				       cluster_size * 2 * M_PI;
			points[i].coords[0] =
				5.0 + 0.5 * cos(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
			points[i].coords[1] =
				2.0 + 0.5 * sin(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
		} else if (i < 3 * cluster_size) {
			/* Cluster 3: centered at (3.5, 5) */
			double angle = (double)(i - 2 * cluster_size) /
				       cluster_size * 2 * M_PI;
			points[i].coords[0] =
				3.5 + 0.5 * cos(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
			points[i].coords[1] =
				5.0 + 0.5 * sin(angle) +
				(rand() / (double)RAND_MAX - 0.5) * 0.2;
		} else {
			/* Random noise points */
			points[i].coords[0] = rand() / (double)RAND_MAX * 7.0;
			points[i].coords[1] = rand() / (double)RAND_MAX * 7.0;
		}
	}
}

/* Print clustering results */
void print_results(cdbscan_point_t *points, int num_points, int num_clusters)
{
	printf("DBSCAN Clustering Results:\n");
	printf("Number of clusters found: %d\n", num_clusters);

	/* Count points in each cluster and noise */
	int *cluster_counts = (int *)calloc(num_clusters + 1, sizeof(int));
	int noise_count = 0;

	for (int i = 0; i < num_points; i++) {
		if (points[i].cluster_id == CDBSCAN_NOISE) {
			noise_count++;
		} else if (points[i].cluster_id >= 0) {
			cluster_counts[points[i].cluster_id]++;
		}
	}

	printf("Noise points: %d\n", noise_count);
	for (int i = 0; i < num_clusters; i++) {
		printf("Cluster %d: %d points\n", i, cluster_counts[i]);
	}

	free(cluster_counts);

	/* Print first 10 points as sample */
	printf("\nSample points (first 10):\n");
	printf("Index\tX\tY\tCluster\n");
	for (int i = 0; i < 10 && i < num_points; i++) {
		printf("%d\t%.2f\t%.2f\t", i, points[i].coords[0],
		       points[i].coords[1]);
		if (points[i].cluster_id == CDBSCAN_NOISE) {
			printf("NOISE\n");
		} else {
			printf("%d\n", points[i].cluster_id);
		}
	}
}

int main()
{
	/* Parameters */
	int num_points = 200;
	int dimensions = 2;
	double eps = 0.5; /* Radius for neighborhood */
	int min_pts = 4; /* Minimum points for core point */

	printf("DBSCAN Clustering Example\n");
	printf("========================\n");
	printf("Number of points: %d\n", num_points);
	printf("Dimensions: %d\n", dimensions);
	printf("Eps (radius): %.2f\n", eps);
	printf("MinPts: %d\n\n", min_pts);

	/* Create points */
	cdbscan_point_t *points = cdbscan_create_points(num_points, dimensions);
	if (!points) {
		fprintf(stderr, "Failed to allocate points\n");
		return 1;
	}

	/* Generate sample data */
	generate_sample_data(points, num_points);

	/* Set up parameters */
	cdbscan_params_t params = { .eps = eps, .min_pts = min_pts };

	/* Run DBSCAN clustering */
	int num_clusters = cdbscan_cluster(points, num_points, params);

	if (num_clusters < 0) {
		fprintf(stderr, "Clustering failed\n");
		cdbscan_free_points(points);
		return 1;
	}

	/* Print results */
	print_results(points, num_points, num_clusters);

	/* Clean up */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);

	return 0;
}

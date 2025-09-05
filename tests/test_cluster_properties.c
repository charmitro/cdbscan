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

/* Test: Cluster Properties - Maximality and Connectivity (Definition 5)
 * 
 * From the paper, a cluster must satisfy:
 * 1) Maximality: if p in C and q is density-reachable from p, then q in C
 * 2) Connectivity: all points in C are density-connected
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cdbscan.h"

void test_cluster_maximality()
{
	printf("Test: Cluster Maximality Property\n");
	printf("==================================\n");

	/* Create scenario where maximality can be tested */
	int num_points = 8;
	cdbscan_point_t *points = cdbscan_create_points(num_points, 2);

	/* Dense group that should form one maximal cluster */
	points[0].coords[0] = 0.0;
	points[0].coords[1] = 0.0; /* Core */
	points[1].coords[0] = 0.2;
	points[1].coords[1] = 0.0; /* Core */
	points[2].coords[0] = 0.1;
	points[2].coords[1] = 0.2; /* Core */
	points[3].coords[0] = 0.4;
	points[3].coords[1] = 0.0; /* Reachable from 1 */
	points[4].coords[0] = 0.3;
	points[4].coords[1] = 0.2; /* Reachable from 2 */
	points[5].coords[0] = -0.2;
	points[5].coords[1] = 0.0; /* Reachable from 0 */

	/* Points just outside reach */
	points[6].coords[0] = 0.7;
	points[6].coords[1] = 0.0; /* Too far */
	points[7].coords[0] = 5.0;
	points[7].coords[1] = 5.0; /* Isolated */

	double eps = 0.25;
	int min_pts = 3;

	printf("Setup: Testing maximality with Eps=%.2f, MinPts=%d\n\n", eps,
	       min_pts);

	/* Run DBSCAN */
	cdbscan_params_t params = { .eps = eps,
				    .min_pts = min_pts,
				    .dist_type = CDBSCAN_DIST_EUCLIDEAN };

	int num_clusters = cdbscan_cluster(points, num_points, params);

	printf("Number of clusters: %d\n\n", num_clusters);

	/* Test maximality: all density-reachable points must be in same cluster */
	int main_cluster = -1;
	for (int i = 0; i <= 5; i++) {
		if (points[i].cluster_id >= 0) {
			if (main_cluster == -1) {
				main_cluster = points[i].cluster_id;
			}
			printf("Point %d: cluster %d ", i,
			       points[i].cluster_id);
			assert(points[i].cluster_id == main_cluster);
			printf("[OK] All reachable points in same cluster\n");
		}
	}

	/* Points 6 and 7 should not be in main cluster */
	printf("\nPoint 6: cluster %2d ", points[6].cluster_id);
	assert(points[6].cluster_id != main_cluster ||
	       points[6].cluster_id == CDBSCAN_NOISE);
	printf("[OK] Not reachable, correctly separated\n");

	printf("Point 7: cluster %2d ", points[7].cluster_id);
	assert(points[7].cluster_id == CDBSCAN_NOISE);
	printf("[OK] Isolated point marked as noise\n");

	printf("\n[PASS] Cluster maximality test passed\n");

	/* Cleanup */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);
}

void test_cluster_connectivity()
{
	printf("\nTest: Cluster Connectivity Property\n");
	printf("====================================\n");

	/* Create two groups that should be separate clusters */
	int num_points = 10;
	cdbscan_point_t *points = cdbscan_create_points(num_points, 2);

	/* Cluster 1 */
	points[0].coords[0] = 0.0;
	points[0].coords[1] = 0.0;
	points[1].coords[0] = 0.2;
	points[1].coords[1] = 0.0;
	points[2].coords[0] = 0.1;
	points[2].coords[1] = 0.2;
	points[3].coords[0] = 0.0;
	points[3].coords[1] = 0.2;
	points[4].coords[0] = -0.1;
	points[4].coords[1] = 0.1;

	/* Cluster 2 - separated from Cluster 1 */
	points[5].coords[0] = 2.0;
	points[5].coords[1] = 0.0;
	points[6].coords[0] = 2.2;
	points[6].coords[1] = 0.0;
	points[7].coords[0] = 2.1;
	points[7].coords[1] = 0.2;
	points[8].coords[0] = 2.0;
	points[8].coords[1] = 0.2;
	points[9].coords[0] = 1.9;
	points[9].coords[1] = 0.1;

	double eps = 0.3;
	int min_pts = 3;

	printf("Setup: Two separate groups with Eps=%.2f, MinPts=%d\n\n", eps,
	       min_pts);

	/* Run DBSCAN */
	cdbscan_params_t params = { .eps = eps,
				    .min_pts = min_pts,
				    .dist_type = CDBSCAN_DIST_EUCLIDEAN };

	int num_clusters = cdbscan_cluster(points, num_points, params);

	printf("Number of clusters: %d\n\n", num_clusters);
	assert(num_clusters == 2);

	/* Verify connectivity within each cluster */
	int cluster1 = points[0].cluster_id;
	int cluster2 = points[5].cluster_id;

	assert(cluster1 != cluster2);
	assert(cluster1 >= 0 && cluster2 >= 0);

	printf("Cluster 1 (points 0-4):\n");
	for (int i = 0; i <= 4; i++) {
		printf("  Point %d: cluster %d ", i, points[i].cluster_id);
		assert(points[i].cluster_id == cluster1);
		printf("[OK] Connected within cluster\n");
	}

	printf("\nCluster 2 (points 5-9):\n");
	for (int i = 5; i <= 9; i++) {
		printf("  Point %d: cluster %d ", i, points[i].cluster_id);
		assert(points[i].cluster_id == cluster2);
		printf("[OK] Connected within cluster\n");
	}

	printf("\n[PASS] Cluster connectivity test passed\n");
	printf("All points within each cluster are density-connected\n");

	/* Cleanup */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);
}

int main()
{
	printf("Testing DBSCAN Cluster Properties (Maximality & Connectivity)\n");
	printf("=============================================================\n\n");

	test_cluster_maximality();
	test_cluster_connectivity();

	printf("\n[SUCCESS] All cluster property tests passed!\n");
	return 0;
}

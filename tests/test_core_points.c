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

/* Test: Core Point Identification (Definition 2 from paper)
 * 
 * From the paper: "A point p is directly density-reachable from a point q 
 * wrt. Eps, MinPts if:
 * 1) p ∈ NEps(q) and
 * 2) |NEps(q)| ≥ MinPts (core point condition)"
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cdbscan.h"

void test_core_point_condition()
{
	printf("Test: Core Point Identification\n");
	printf("================================\n");

	/* Create a simple test case with known core and non-core points */
	int num_points = 10;
	cdbscan_point_t *points = cdbscan_create_points(num_points, 2);

	/* Cluster of 5 points at origin (should be core points) */
	points[0].coords[0] = 0.0;
	points[0].coords[1] = 0.0;
	points[1].coords[0] = 0.1;
	points[1].coords[1] = 0.0;
	points[2].coords[0] = 0.0;
	points[2].coords[1] = 0.1;
	points[3].coords[0] = -0.1;
	points[3].coords[1] = 0.0;
	points[4].coords[0] = 0.0;
	points[4].coords[1] = -0.1;

	/* Cluster of 3 points at (5,5) (not enough for core with MinPts=4) */
	points[5].coords[0] = 5.0;
	points[5].coords[1] = 5.0;
	points[6].coords[0] = 5.1;
	points[6].coords[1] = 5.0;
	points[7].coords[0] = 5.0;
	points[7].coords[1] = 5.1;

	/* Two isolated points (definitely not core) */
	points[8].coords[0] = 10.0;
	points[8].coords[1] = 10.0;
	points[9].coords[0] = -10.0;
	points[9].coords[1] = -10.0;

	/* Test parameters */
	double eps = 0.3;
	int min_pts = 4;

	printf("Testing with Eps=%.2f, MinPts=%d\n\n", eps, min_pts);

	/* Check neighborhoods for each point */
	int *neighbors = (int *)malloc(num_points * sizeof(int));

	for (int i = 0; i < num_points; i++) {
		int neighbor_count = cdbscan_region_query(points, num_points, i,
							  eps, neighbors);

		printf("Point %d (%.1f, %.1f): %d neighbors -> ", i,
		       points[i].coords[0], points[i].coords[1],
		       neighbor_count);

		if (neighbor_count >= min_pts) {
			printf("CORE POINT\n");

			/* Verify: Points 0-4 should be core points */
			if (i <= 4) {
				printf("  [OK] Correctly identified as core\n");
			} else {
				printf("  [ERROR] Should not be core\n");
				assert(0);
			}
		} else {
			printf("NOT CORE\n");

			/* Verify: Points 5-9 should not be core points */
			if (i >= 5) {
				printf("  [OK] Correctly identified as non-core\n");
			} else {
				printf("  [ERROR] Should be core\n");
				assert(0);
			}
		}
	}

	/* Run DBSCAN and verify clustering */
	cdbscan_params_t params = { .eps = eps,
				    .min_pts = min_pts,
				    .dist_type = CDBSCAN_DIST_EUCLIDEAN };

	int num_clusters = cdbscan_cluster(points, num_points, params);

	printf("\n=== DBSCAN Results ===\n");
	printf("Number of clusters: %d\n", num_clusters);

	/* Verify that core points (0-4) are in the same cluster */
	int cluster_of_core = points[0].cluster_id;
	assert(cluster_of_core >= 0); /* Should not be noise */

	for (int i = 1; i <= 4; i++) {
		assert(points[i].cluster_id == cluster_of_core);
		printf("Point %d: cluster %d [OK]\n", i, points[i].cluster_id);
	}

	printf("\n[PASS] Core point test PASSED\n");

	/* Cleanup */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);
	free(neighbors);
}

int main()
{
	printf("Testing DBSCAN Core Point Specification\n");
	printf("========================================\n\n");

	test_core_point_condition();

	printf("\n[SUCCESS] All core point tests passed!\n");
	return 0;
}

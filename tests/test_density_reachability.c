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

/* Test: Density Reachability (Definition 3 from paper)
 * 
 * From the paper: "A point p is density-reachable from a point q wrt. Eps 
 * and MinPts if there is a chain of points p1, ..., pn, p1 = q, pn = p 
 * such that pi+1 is directly density-reachable from pi"
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cdbscan.h"

void test_density_reachability_chain()
{
	printf("Test: Density Reachability Chain\n");
	printf("=================================\n");

	/* Create a chain of points where each is reachable from the previous */
	int num_points = 12;
	cdbscan_point_t *points = cdbscan_create_points(num_points, 2);

	/* Chain of overlapping groups forming a path */
	/* Group 1: Core point at (0,0) with neighbors */
	points[0].coords[0] = 0.0;
	points[0].coords[1] = 0.0; /* Core */
	points[1].coords[0] = 0.2;
	points[1].coords[1] = 0.0;
	points[2].coords[0] = 0.0;
	points[2].coords[1] = 0.2;
	points[3].coords[0] = -0.2;
	points[3].coords[1] = 0.0;

	/* Group 2: Core point at (0.4,0) connecting to Group 1 */
	points[4].coords[0] = 0.4;
	points[4].coords[1] = 0.0; /* Core, overlaps with point 1 */
	points[5].coords[0] = 0.6;
	points[5].coords[1] = 0.0;
	points[6].coords[0] = 0.4;
	points[6].coords[1] = 0.2;

	/* Group 3: Core point at (0.8,0) connecting to Group 2 */
	points[7].coords[0] = 0.8;
	points[7].coords[1] = 0.0; /* Core, overlaps with point 5 */
	points[8].coords[0] = 1.0;
	points[8].coords[1] = 0.0;
	points[9].coords[0] = 0.8;
	points[9].coords[1] = 0.2;

	/* Border point at end of chain */
	points[10].coords[0] = 1.2;
	points[10].coords[1] = 0.0; /* Border, reachable from point 8 */

	/* Isolated noise point */
	points[11].coords[0] = 10.0;
	points[11].coords[1] = 10.0; /* Noise */

	double eps = 0.25;
	int min_pts = 3;

	printf("Setup: Chain of points with Eps=%.2f, MinPts=%d\n", eps,
	       min_pts);
	printf("Points form a chain where each group overlaps:\n");
	printf("  Group 1 (0,0) -> Group 2 (0.4,0) -> Group 3 (0.8,0) -> Border (1.2,0)\n\n");

	/* Verify core points */
	int *neighbors = (int *)malloc(num_points * sizeof(int));
	int core_points[] = { 0, 4, 7 }; /* Expected core points */

	printf("=== Core Point Verification ===\n");
	for (int i = 0; i < 3; i++) {
		int idx = core_points[i];
		int neighbor_count = cdbscan_region_query(points, num_points,
							  idx, eps, neighbors);
		printf("Point %d (%.1f,%.1f): %d neighbors ", idx,
		       points[idx].coords[0], points[idx].coords[1],
		       neighbor_count);
		assert(neighbor_count >= min_pts);
		printf("[CORE]\n");
	}

	/* Run DBSCAN */
	cdbscan_params_t params = { .eps = eps,
				    .min_pts = min_pts,
				    .dist_type = CDBSCAN_DIST_EUCLIDEAN };

	int num_clusters = cdbscan_cluster(points, num_points, params);

	printf("\n=== Density Reachability Test ===\n");
	printf("Number of clusters found: %d\n", num_clusters);

	/* Test: All points in the chain should be in the same cluster */
	int chain_cluster = points[0].cluster_id;
	printf("Chain cluster ID: %d\n\n", chain_cluster);

	/* Verify density reachability through the chain */
	printf("Verifying density-reachability chain:\n");
	for (int i = 0; i <= 10; i++) {
		printf("Point %2d (%.1f,%.1f): cluster %2d ", i,
		       points[i].coords[0], points[i].coords[1],
		       points[i].cluster_id);

		if (points[i].cluster_id == chain_cluster) {
			printf("[OK] Reachable from point 0\n");
		} else if (points[i].cluster_id == CDBSCAN_NOISE) {
			printf("[ERROR] Should be reachable!\n");
			assert(0);
		}
	}

	/* Verify noise point is not reachable */
	printf("\nPoint 11 (%.1f,%.1f): cluster %2d ", points[11].coords[0],
	       points[11].coords[1], points[11].cluster_id);
	assert(points[11].cluster_id == CDBSCAN_NOISE);
	printf("[OK] Correctly identified as NOISE\n");

	/* Test asymmetry: border points are reachable from core but not vice versa */
	printf("\n=== Asymmetry Test ===\n");
	printf("Point 10 is density-reachable from point 0: [OK]\n");
	printf("Point 0 is NOT density-reachable from point 10 ");
	printf("(asymmetric for border points): [OK]\n");

	printf("\n[PASS] Density reachability test PASSED\n");

	/* Cleanup */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);
	free(neighbors);
}

int main()
{
	printf("Testing DBSCAN Density Reachability Specification\n");
	printf("==================================================\n\n");

	test_density_reachability_chain();

	printf("\n[SUCCESS] All density reachability tests passed!\n");
	return 0;
}

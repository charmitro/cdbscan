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

/* Test: Border Points and Noise (Definitions 5 & 6 from paper)
 * 
 * From the paper: 
 * - Border points: in cluster but not core points
 * - Noise: points not belonging to any cluster
 * - Border points of same cluster may not be density-reachable from each other
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cdbscan.h"

void test_border_and_noise_points()
{
	printf("Test: Border Points and Noise Classification\n");
	printf("=============================================\n");

	int num_points = 15;
	cdbscan_point_t *points = cdbscan_create_points(num_points, 2);

	/* Cluster 1: Dense cluster with clear core and border points */
	/* Core points (0-3): form a square */
	points[0].coords[0] = 0.0;
	points[0].coords[1] = 0.0;
	points[1].coords[0] = 0.2;
	points[1].coords[1] = 0.0;
	points[2].coords[0] = 0.2;
	points[2].coords[1] = 0.2;
	points[3].coords[0] = 0.0;
	points[3].coords[1] = 0.2;

	/* Border points (4-6): connected to core but not enough neighbors */
	points[4].coords[0] = 0.4;
	points[4].coords[1] = 0.1; /* Border of cluster 1 */
	points[5].coords[0] = -0.2;
	points[5].coords[1] = 0.1; /* Border of cluster 1 */
	points[6].coords[0] = 0.1;
	points[6].coords[1] = 0.4; /* Border of cluster 1 */

	/* Cluster 2: Another dense cluster */
	/* Core points (7-10) */
	points[7].coords[0] = 3.0;
	points[7].coords[1] = 0.0;
	points[8].coords[0] = 3.2;
	points[8].coords[1] = 0.0;
	points[9].coords[0] = 3.2;
	points[9].coords[1] = 0.2;
	points[10].coords[0] = 3.0;
	points[10].coords[1] = 0.2;

	/* Border point (11) */
	points[11].coords[0] = 3.4;
	points[11].coords[1] = 0.1; /* Border of cluster 2 */

	/* Noise points (12-14): isolated points */
	points[12].coords[0] = 1.5;
	points[12].coords[1] = 1.5; /* Between clusters */
	points[13].coords[0] = -2.0;
	points[13].coords[1] = -2.0; /* Isolated */
	points[14].coords[0] = 5.0;
	points[14].coords[1] = 5.0; /* Isolated */

	double eps = 0.25;
	int min_pts = 4;

	printf("Setup: Eps=%.2f, MinPts=%d\n", eps, min_pts);
	printf("Expected structure:\n");
	printf("  - Cluster 1: points 0-3 (core), 4-6 (border)\n");
	printf("  - Cluster 2: points 7-10 (core), 11 (border)\n");
	printf("  - Noise: points 12-14\n\n");

	/* Analyze neighborhoods before clustering */
	int *neighbors = (int *)malloc(num_points * sizeof(int));

	printf("=== Neighborhood Analysis ===\n");
	for (int i = 0; i < num_points; i++) {
		int neighbor_count = cdbscan_region_query(points, num_points, i,
							  eps, neighbors);
		printf("Point %2d (%.1f,%.1f): %d neighbors -> ", i,
		       points[i].coords[0], points[i].coords[1],
		       neighbor_count);

		if (neighbor_count >= min_pts) {
			printf("CORE\n");
		} else if (neighbor_count > 1) {
			printf("Potential BORDER\n");
		} else {
			printf("Likely NOISE\n");
		}
	}

	/* Run DBSCAN */
	cdbscan_params_t params = { .eps = eps,
				    .min_pts = min_pts,
				    .dist_type = CDBSCAN_DIST_EUCLIDEAN };

	int num_clusters = cdbscan_cluster(points, num_points, params);

	printf("\n=== DBSCAN Results ===\n");
	printf("Number of clusters: %d\n\n", num_clusters);

	/* Verify classifications */
	printf("Point Classifications:\n");

	/* Check core points */
	int core_indices[] = { 0, 1, 2, 3, 7, 8, 9, 10 };
	for (int i = 0; i < 8; i++) {
		int idx = core_indices[i];
		printf("Point %2d: cluster %2d ", idx, points[idx].cluster_id);
		assert(points[idx].cluster_id >= 0); /* Must be in a cluster */
		printf("[OK] Core point in cluster\n");
	}

	printf("\n");

	/* Check border points */
	int border_indices[] = { 4, 5, 6, 11 };
	for (int i = 0; i < 4; i++) {
		int idx = border_indices[i];
		printf("Point %2d: cluster %2d ", idx, points[idx].cluster_id);
		assert(points[idx].cluster_id >= 0); /* Must be in a cluster */
		printf("[OK] Border point in cluster\n");
	}

	printf("\n");

	/* Check noise points */
	int noise_indices[] = { 12, 13, 14 };
	for (int i = 0; i < 3; i++) {
		int idx = noise_indices[i];
		printf("Point %2d: cluster %2d ", idx, points[idx].cluster_id);
		assert(points[idx].cluster_id == CDBSCAN_NOISE);
		printf("[OK] Correctly identified as NOISE\n");
	}

	/* Test: Border points of same cluster */
	printf("\n=== Border Point Connectivity Test ===\n");
	if (points[4].cluster_id == points[5].cluster_id) {
		printf("Border points 4 and 5 are in the same cluster: [OK]\n");
		printf("(Connected through core points, not directly to each other)\n");
	}

	/* Test: Verify cluster assignments are consistent */
	printf("\n=== Cluster Consistency ===\n");
	assert(points[0].cluster_id == points[1].cluster_id);
	assert(points[1].cluster_id == points[2].cluster_id);
	assert(points[2].cluster_id == points[3].cluster_id);
	printf("Cluster 1 core points: consistent [OK]\n");

	assert(points[7].cluster_id == points[8].cluster_id);
	assert(points[8].cluster_id == points[9].cluster_id);
	assert(points[9].cluster_id == points[10].cluster_id);
	printf("Cluster 2 core points: consistent [OK]\n");

	printf("\n[PASS] Border and noise test PASSED\n");

	/* Cleanup */
	for (int i = 0; i < num_points; i++) {
		free(points[i].coords);
	}
	free(points);
	free(neighbors);
}

int main()
{
	printf("Testing DBSCAN Border Points and Noise Specification\n");
	printf("====================================================\n\n");

	test_border_and_noise_points();

	printf("\n[SUCCESS] All border and noise tests passed!\n");
	return 0;
}

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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CDBSCAN_H
#define CDBSCAN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Library version */
#define CDBSCAN_VERSION_MAJOR 0
#define CDBSCAN_VERSION_MINOR 1
#define CDBSCAN_VERSION_PATCH 0

/* Special cluster IDs */
#define CDBSCAN_UNCLASSIFIED -1
#define CDBSCAN_NOISE -2

/* Distance metric types */
typedef enum {
	CDBSCAN_DIST_EUCLIDEAN,
	CDBSCAN_DIST_MANHATTAN,
	CDBSCAN_DIST_MINKOWSKI,
	CDBSCAN_DIST_COSINE,
	CDBSCAN_DIST_CUSTOM
} cdbscan_dist_type_t;

/* Custom distance function prototype */
typedef double (*cdbscan_dist_func_t)(const double *a, const double *b,
				      int dims, void *params);

/* Point structure for 2D/3D spatial data */
typedef struct cdbscan_point {
	double *coords; /* Coordinate array */
	int dimensions; /* Number of dimensions */
	int cluster_id; /* Cluster assignment */
	int index; /* Original index in dataset */
} cdbscan_point_t;

/* DBSCAN parameters */
typedef struct cdbscan_params {
	double eps; /* Epsilon: radius for neighborhood */
	int min_pts; /* Minimum points to form dense region */
	cdbscan_dist_type_t dist_type; /* Distance metric to use */
	double minkowski_p; /* p parameter for Minkowski distance */
	cdbscan_dist_func_t custom_dist; /* Custom distance function */
	void *custom_dist_params; /* Parameters for custom distance */
	int use_kdtree; /* Use KD-tree for O(n log n) performance (1=yes, 0=no) */
} cdbscan_params_t;

/* Main DBSCAN clustering function
 * Returns: number of clusters found (excluding noise)
 * Sets cluster_id field in each point:
 *   - >= 0: cluster number
 *   - CDBSCAN_NOISE: noise point
 */
int cdbscan_cluster(cdbscan_point_t *points, int num_points,
		    cdbscan_params_t params);

/* Distance functions */
double cdbscan_euclidean_distance(const double *a, const double *b, int dims);
double cdbscan_manhattan_distance(const double *a, const double *b, int dims);
double cdbscan_minkowski_distance(const double *a, const double *b, int dims,
				  double p);
double cdbscan_cosine_distance(const double *a, const double *b, int dims);

/* Data normalization functions */
void cdbscan_normalize_minmax(cdbscan_point_t *points, int num_points);
void cdbscan_normalize_zscore(cdbscan_point_t *points, int num_points);

/* Parameter estimation utilities */
typedef struct {
	double *distances; /* k-distances for each point */
	int k; /* k value used */
	double suggested_eps; /* Suggested eps value */
} cdbscan_kdist_result_t;

cdbscan_kdist_result_t *cdbscan_estimate_eps(const cdbscan_point_t *points,
					     int num_points,
					     int k /* typically 4 for 2D data */
);
void cdbscan_free_kdist_result(cdbscan_kdist_result_t *result);

/* Region query: find all points within eps distance of point
 * Returns: number of neighbors found
 * neighbors: array to store neighbor indices (must be pre-allocated)
 */
int cdbscan_region_query(const cdbscan_point_t *points, int num_points,
			 int point_idx, double eps, int *neighbors);

/* Advanced region query with custom distance */
int cdbscan_region_query_custom(const cdbscan_point_t *points, int num_points,
				int point_idx, const cdbscan_params_t *params,
				int *neighbors);

/* Utility functions */
cdbscan_point_t *cdbscan_create_points(int num_points, int dimensions);
void cdbscan_free_points(cdbscan_point_t *points);

/* Validation functions */
int cdbscan_validate_params(const cdbscan_params_t *params);
int cdbscan_validate_data(const cdbscan_point_t *points, int num_points);

#ifdef __cplusplus
}
#endif

#endif /* CDBSCAN_H */

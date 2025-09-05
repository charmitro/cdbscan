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

#include "cdbscan.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Internal comparison function for qsort */
static int compare_doubles(const void *a, const void *b)
{
	double diff = *(double *)a - *(double *)b;
	return (diff > 0) - (diff < 0);
}

/* Distance metric implementations */
double cdbscan_euclidean_distance(const double *a, const double *b, int dims)
{
	if (!a || !b || dims <= 0)
		return -1.0;

	double sum = 0.0;
	for (int i = 0; i < dims; i++) {
		double diff = a[i] - b[i];
		sum += diff * diff;
	}
	return sqrt(sum);
}

double cdbscan_manhattan_distance(const double *a, const double *b, int dims)
{
	if (!a || !b || dims <= 0)
		return -1.0;

	double sum = 0.0;
	for (int i = 0; i < dims; i++) {
		sum += fabs(a[i] - b[i]);
	}
	return sum;
}

double cdbscan_minkowski_distance(const double *a, const double *b, int dims,
				  double p)
{
	if (!a || !b || dims <= 0 || p <= 0)
		return -1.0;

	double sum = 0.0;
	for (int i = 0; i < dims; i++) {
		sum += pow(fabs(a[i] - b[i]), p);
	}
	return pow(sum, 1.0 / p);
}

double cdbscan_cosine_distance(const double *a, const double *b, int dims)
{
	if (!a || !b || dims <= 0)
		return -1.0;

	double dot = 0.0, norm_a = 0.0, norm_b = 0.0;
	for (int i = 0; i < dims; i++) {
		dot += a[i] * b[i];
		norm_a += a[i] * a[i];
		norm_b += b[i] * b[i];
	}

	if (norm_a == 0.0 || norm_b == 0.0)
		return 2.0; /* Maximum distance */

	double similarity = dot / (sqrt(norm_a) * sqrt(norm_b));
	return 1.0 - similarity; /* Convert similarity to distance */
}

/* Internal distance calculation based on params */
static double calculate_distance(const double *a, const double *b, int dims,
				 const cdbscan_params_t *params)
{
	switch (params->dist_type) {
	case CDBSCAN_DIST_EUCLIDEAN:
		return cdbscan_euclidean_distance(a, b, dims);
	case CDBSCAN_DIST_MANHATTAN:
		return cdbscan_manhattan_distance(a, b, dims);
	case CDBSCAN_DIST_MINKOWSKI:
		return cdbscan_minkowski_distance(a, b, dims,
						  params->minkowski_p);
	case CDBSCAN_DIST_COSINE:
		return cdbscan_cosine_distance(a, b, dims);
	case CDBSCAN_DIST_CUSTOM:
		if (params->custom_dist) {
			return params->custom_dist(a, b, dims,
						   params->custom_dist_params);
		}
		break;
	}
	return -1.0;
}

/* KD-tree implementation for O(n log n) performance */
typedef struct kdtree_node {
	int point_idx; /* Index of point in original array */
	struct kdtree_node *left;
	struct kdtree_node *right;
	int split_dim; /* Dimension used for splitting at this node */
} kdtree_node_t;

typedef struct {
	kdtree_node_t *root;
	const cdbscan_point_t *points; /* Reference to original points */
	int num_points;
	int dimensions;
} kdtree_t;

/* Helper: Select median and partition array */
static int partition(int *indices, const cdbscan_point_t *points, int left,
		     int right, int dim)
{
	int pivot_idx = (left + right) / 2;
	int pivot_point = indices[pivot_idx];
	double pivot_val = points[pivot_point].coords[dim];

	/* Move pivot to end */
	int temp = indices[pivot_idx];
	indices[pivot_idx] = indices[right];
	indices[right] = temp;

	int store_idx = left;
	for (int i = left; i < right; i++) {
		if (points[indices[i]].coords[dim] < pivot_val) {
			temp = indices[store_idx];
			indices[store_idx] = indices[i];
			indices[i] = temp;
			store_idx++;
		}
	}

	/* Move pivot to its final place */
	temp = indices[store_idx];
	indices[store_idx] = indices[right];
	indices[right] = temp;

	return store_idx;
}

/* Helper: Perform nth_element partitioning (like C++ std::nth_element) */
static void nth_element(int *indices, const cdbscan_point_t *points, int left,
			int right, int n, int dim)
{
	while (left < right) {
		int pivot_idx = partition(indices, points, left, right, dim);

		if (pivot_idx == n) {
			return;
		} else if (pivot_idx > n) {
			right = pivot_idx - 1;
		} else {
			left = pivot_idx + 1;
		}
	}
}

/* Build KD-tree recursively */
static kdtree_node_t *kdtree_build_recursive(int *indices, int num_indices,
					     const cdbscan_point_t *points,
					     int depth, int dimensions)
{
	if (num_indices <= 0)
		return NULL;

	kdtree_node_t *node = (kdtree_node_t *)calloc(1, sizeof(kdtree_node_t));
	if (!node)
		return NULL;

	if (num_indices == 1) {
		node->point_idx = indices[0];
		node->split_dim = depth % dimensions;
		return node;
	}

	/* Choose splitting dimension (cycle through dimensions) */
	int split_dim = depth % dimensions;
	node->split_dim = split_dim;

	/* Find median position */
	int median_idx = num_indices / 2;

	/* Partition array so median is at correct position */
	nth_element(indices, points, 0, num_indices - 1, median_idx, split_dim);

	node->point_idx = indices[median_idx];

	/* Recursively build left and right subtrees */
	node->left = kdtree_build_recursive(indices, median_idx, points,
					    depth + 1, dimensions);
	node->right = kdtree_build_recursive(indices + median_idx + 1,
					     num_indices - median_idx - 1,
					     points, depth + 1, dimensions);

	return node;
}

/* Build KD-tree from points */
static kdtree_t *kdtree_build(const cdbscan_point_t *points, int num_points)
{
	if (!points || num_points <= 0)
		return NULL;

	kdtree_t *tree = (kdtree_t *)calloc(1, sizeof(kdtree_t));
	if (!tree)
		return NULL;

	/* Create array of indices */
	int *indices = (int *)malloc(num_points * sizeof(int));
	if (!indices) {
		free(tree);
		return NULL;
	}

	for (int i = 0; i < num_points; i++) {
		indices[i] = i;
	}

	tree->points = points;
	tree->num_points = num_points;
	tree->dimensions = points[0].dimensions;
	tree->root = kdtree_build_recursive(indices, num_points, points, 0,
					    tree->dimensions);

	free(indices);
	return tree;
}

/* Free KD-tree recursively */
static void kdtree_free_recursive(kdtree_node_t *node)
{
	if (!node)
		return;
	kdtree_free_recursive(node->left);
	kdtree_free_recursive(node->right);
	free(node);
}

/* Free KD-tree */
static void kdtree_free(kdtree_t *tree)
{
	if (!tree)
		return;
	kdtree_free_recursive(tree->root);
	free(tree);
}

/* Range query: find all points within eps distance */
static void kdtree_range_query_recursive(const kdtree_node_t *node,
					 const cdbscan_point_t *query_point,
					 double eps, double eps_squared,
					 const cdbscan_point_t *points,
					 int *neighbors, int *count,
					 int dimensions)
{
	if (!node)
		return;

	const cdbscan_point_t *node_point = &points[node->point_idx];

	/* Calculate actual Euclidean distance */
	double dist = cdbscan_euclidean_distance(
		query_point->coords, node_point->coords, dimensions);

	/* If within range, add to neighbors */
	if (dist <= eps) {
		neighbors[(*count)++] = node->point_idx;
	}

	/* Get splitting dimension and value */
	int split_dim = node->split_dim;
	double split_val = node_point->coords[split_dim];
	double query_val = query_point->coords[split_dim];
	double diff = query_val - split_val;

	/* Determine which subtree to search first */
	kdtree_node_t *first_child = (diff < 0) ? node->left : node->right;
	kdtree_node_t *second_child = (diff < 0) ? node->right : node->left;

	/* Search the closer subtree first */
	kdtree_range_query_recursive(first_child, query_point, eps, eps_squared,
				     points, neighbors, count, dimensions);

	/* Only search the other subtree if it could contain points within eps */
	if (fabs(diff) <= eps) {
		kdtree_range_query_recursive(second_child, query_point, eps,
					     eps_squared, points, neighbors,
					     count, dimensions);
	}
}

/* Helper: Compare function for sorting integers */
static int compare_ints(const void *a, const void *b)
{
	return *(int *)a - *(int *)b;
}

/* KD-tree range query */
static int kdtree_range_query(const kdtree_t *tree, int query_idx, double eps,
			      int *neighbors)
{
	if (!tree || !tree->root || !neighbors)
		return 0;

	int count = 0;
	double eps_squared = eps * eps;
	const cdbscan_point_t *query_point = &tree->points[query_idx];

	kdtree_range_query_recursive(tree->root, query_point, eps, eps_squared,
				     tree->points, neighbors, &count,
				     tree->dimensions);

	/* Sort neighbors to ensure consistent ordering */
	if (count > 0) {
		qsort(neighbors, count, sizeof(int), compare_ints);
	}

	return count;
}

/* Data normalization functions */
void cdbscan_normalize_minmax(cdbscan_point_t *points, int num_points)
{
	if (!points || num_points <= 0)
		return;

	int dims = points[0].dimensions;
	double *min_vals = (double *)calloc(dims, sizeof(double));
	double *max_vals = (double *)calloc(dims, sizeof(double));

	if (!min_vals || !max_vals) {
		free(min_vals);
		free(max_vals);
		return;
	}

	/* Initialize min/max */
	for (int d = 0; d < dims; d++) {
		min_vals[d] = DBL_MAX;
		max_vals[d] = -DBL_MAX;
	}

	/* Find min/max for each dimension */
	for (int i = 0; i < num_points; i++) {
		for (int d = 0; d < dims; d++) {
			if (points[i].coords[d] < min_vals[d]) {
				min_vals[d] = points[i].coords[d];
			}
			if (points[i].coords[d] > max_vals[d]) {
				max_vals[d] = points[i].coords[d];
			}
		}
	}

	/* Normalize to [0, 1] */
	for (int i = 0; i < num_points; i++) {
		for (int d = 0; d < dims; d++) {
			double range = max_vals[d] - min_vals[d];
			if (range > 0) {
				points[i].coords[d] =
					(points[i].coords[d] - min_vals[d]) /
					range;
			} else {
				points[i].coords[d] = 0.0;
			}
		}
	}

	free(min_vals);
	free(max_vals);
}

void cdbscan_normalize_zscore(cdbscan_point_t *points, int num_points)
{
	if (!points || num_points <= 0)
		return;

	int dims = points[0].dimensions;
	double *means = (double *)calloc(dims, sizeof(double));
	double *stdevs = (double *)calloc(dims, sizeof(double));

	if (!means || !stdevs) {
		free(means);
		free(stdevs);
		return;
	}

	/* Calculate means */
	for (int i = 0; i < num_points; i++) {
		for (int d = 0; d < dims; d++) {
			means[d] += points[i].coords[d];
		}
	}
	for (int d = 0; d < dims; d++) {
		means[d] /= num_points;
	}

	/* Calculate standard deviations */
	for (int i = 0; i < num_points; i++) {
		for (int d = 0; d < dims; d++) {
			double diff = points[i].coords[d] - means[d];
			stdevs[d] += diff * diff;
		}
	}
	for (int d = 0; d < dims; d++) {
		stdevs[d] = sqrt(stdevs[d] / num_points);
	}

	/* Normalize using z-score */
	for (int i = 0; i < num_points; i++) {
		for (int d = 0; d < dims; d++) {
			if (stdevs[d] > 0) {
				points[i].coords[d] =
					(points[i].coords[d] - means[d]) /
					stdevs[d];
			} else {
				points[i].coords[d] = 0.0;
			}
		}
	}

	free(means);
	free(stdevs);
}

/* Parameter estimation - k-dist graph for eps selection */
cdbscan_kdist_result_t *cdbscan_estimate_eps(const cdbscan_point_t *points,
					     int num_points, int k)
{
	if (!points || num_points <= 0 || k <= 0 || k >= num_points) {
		return NULL;
	}

	cdbscan_kdist_result_t *result = (cdbscan_kdist_result_t *)malloc(
		sizeof(cdbscan_kdist_result_t));
	if (!result)
		return NULL;

	result->distances = (double *)malloc(num_points * sizeof(double));
	if (!result->distances) {
		free(result);
		return NULL;
	}

	result->k = k;

	/* For each point, find k-th nearest neighbor distance */
	double *temp_dists = (double *)malloc(num_points * sizeof(double));
	if (!temp_dists) {
		free(result->distances);
		free(result);
		return NULL;
	}

	for (int i = 0; i < num_points; i++) {
		/* Calculate distances to all other points */
		int dist_count = 0;
		for (int j = 0; j < num_points; j++) {
			if (i != j) {
				temp_dists[dist_count++] =
					cdbscan_euclidean_distance(
						points[i].coords,
						points[j].coords,
						points[i].dimensions);
			}
		}

		/* Sort distances */
		qsort(temp_dists, dist_count, sizeof(double), compare_doubles);

		/* Store k-th distance (k-1 because array is 0-indexed) */
		result->distances[i] = temp_dists[k - 1];
	}

	/* Sort k-distances in descending order for graph */
	memcpy(temp_dists, result->distances, num_points * sizeof(double));
	qsort(temp_dists, num_points, sizeof(double), compare_doubles);

	/* Find the "elbow" - simplified: use 95th percentile */
	int elbow_idx = (int)(0.95 * num_points);
	result->suggested_eps = temp_dists[elbow_idx];

	free(temp_dists);

	return result;
}

void cdbscan_free_kdist_result(cdbscan_kdist_result_t *result)
{
	if (result) {
		free(result->distances);
		free(result);
	}
}

/* Validation functions */
int cdbscan_validate_params(const cdbscan_params_t *params)
{
	if (!params)
		return 0;
	if (params->eps <= 0)
		return 0;
	if (params->min_pts <= 0)
		return 0;

	if (params->dist_type == CDBSCAN_DIST_MINKOWSKI &&
	    params->minkowski_p <= 0) {
		return 0;
	}

	if (params->dist_type == CDBSCAN_DIST_CUSTOM && !params->custom_dist) {
		return 0;
	}

	return 1;
}

int cdbscan_validate_data(const cdbscan_point_t *points, int num_points)
{
	if (!points || num_points <= 0)
		return 0;

	int dims = points[0].dimensions;
	if (dims <= 0)
		return 0;

	for (int i = 0; i < num_points; i++) {
		if (!points[i].coords)
			return 0;
		if (points[i].dimensions != dims)
			return 0;

		/* Check for NaN or Inf */
		for (int d = 0; d < dims; d++) {
			if (isnan(points[i].coords[d]) ||
			    isinf(points[i].coords[d])) {
				return 0;
			}
		}
	}

	return 1;
}

/* Region query implementations */
int cdbscan_region_query(const cdbscan_point_t *points, int num_points,
			 int point_idx, double eps, int *neighbors)
{
	if (!points || !neighbors || point_idx < 0 || point_idx >= num_points) {
		return 0;
	}

	int neighbor_count = 0;
	const cdbscan_point_t *query_point = &points[point_idx];

	for (int i = 0; i < num_points; i++) {
		double dist = cdbscan_euclidean_distance(
			query_point->coords, points[i].coords,
			query_point->dimensions);
		if (dist >= 0 && dist <= eps) {
			neighbors[neighbor_count++] = i;
		}
	}

	return neighbor_count;
}

int cdbscan_region_query_custom(const cdbscan_point_t *points, int num_points,
				int point_idx, const cdbscan_params_t *params,
				int *neighbors)
{
	if (!points || !neighbors || !params || point_idx < 0 ||
	    point_idx >= num_points) {
		return 0;
	}

	int neighbor_count = 0;
	const cdbscan_point_t *query_point = &points[point_idx];

	for (int i = 0; i < num_points; i++) {
		double dist = calculate_distance(query_point->coords,
						 points[i].coords,
						 query_point->dimensions,
						 params);
		if (dist >= 0 && dist <= params->eps) {
			neighbors[neighbor_count++] = i;
		}
	}

	return neighbor_count;
}

/* Forward declaration for internal function */
static int expand_cluster(cdbscan_point_t *points, int num_points,
			  int point_idx, int cluster_id,
			  const cdbscan_params_t *params, int *neighbors,
			  int *seeds, int *seed_size);

/* Forward declaration for KD-tree version */
static int expand_cluster_kdtree(cdbscan_point_t *points, int num_points,
				 int point_idx, int cluster_id,
				 const cdbscan_params_t *params,
				 const kdtree_t *tree, int *neighbors,
				 int *seeds, int *seed_size);

/* Main DBSCAN clustering algorithm */
int cdbscan_cluster(cdbscan_point_t *points, int num_points,
		    cdbscan_params_t params)
{
	/* Validate inputs */
	if (!cdbscan_validate_params(&params))
		return -1;
	if (!cdbscan_validate_data(points, num_points))
		return -1;

	/* Initialize all points as UNCLASSIFIED */
	for (int i = 0; i < num_points; i++) {
		points[i].cluster_id = CDBSCAN_UNCLASSIFIED;
		points[i].index = i;
	}

	/* Allocate working arrays */
	int *neighbors = (int *)malloc(num_points * sizeof(int));
	int *seeds = (int *)malloc(num_points * sizeof(int));
	if (!neighbors || !seeds) {
		free(neighbors);
		free(seeds);
		return -1;
	}

	/* Build KD-tree if requested and using Euclidean distance */
	kdtree_t *tree = NULL;
	if (params.use_kdtree && params.dist_type == CDBSCAN_DIST_EUCLIDEAN) {
		tree = kdtree_build(points, num_points);
		if (!tree) {
			/* Fall back to brute force if tree building fails */
			params.use_kdtree = 0;
		}
	}

	int cluster_id = 0;

	/* Process each point */
	for (int i = 0; i < num_points; i++) {
		if (points[i].cluster_id != CDBSCAN_UNCLASSIFIED) {
			continue; /* Already processed */
		}

		/* Find neighbors using KD-tree or brute force */
		int neighbor_count;
		if (tree) {
			neighbor_count = kdtree_range_query(tree, i, params.eps,
							    neighbors);
		} else {
			neighbor_count = cdbscan_region_query_custom(
				points, num_points, i, &params, neighbors);
		}

		if (neighbor_count < params.min_pts) {
			/* Mark as noise (may be changed later if it's a border point) */
			points[i].cluster_id = CDBSCAN_NOISE;
		} else {
			/* Core point - start a new cluster */
			int seed_size = 0;
			if (tree) {
				if (expand_cluster_kdtree(points, num_points, i,
							  cluster_id, &params,
							  tree, neighbors,
							  seeds, &seed_size)) {
					cluster_id++;
				}
			} else {
				if (expand_cluster(points, num_points, i,
						   cluster_id, &params,
						   neighbors, seeds,
						   &seed_size)) {
					cluster_id++;
				}
			}
		}
	}

	/* Clean up */
	if (tree) {
		kdtree_free(tree);
	}
	free(neighbors);
	free(seeds);

	return cluster_id; /* Return number of clusters found */
}

/* Expand cluster from a core point */
static int expand_cluster(cdbscan_point_t *points, int num_points,
			  int point_idx, int cluster_id,
			  const cdbscan_params_t *params, int *neighbors,
			  int *seeds, int *seed_size)
{
	/* Get initial seeds from region query */
	*seed_size = cdbscan_region_query_custom(points, num_points, point_idx,
						 params, seeds);

	if (*seed_size < params->min_pts) {
		/* Not a core point */
		points[point_idx].cluster_id = CDBSCAN_NOISE;
		return 0;
	}

	/* Assign cluster ID to all points in seeds */
	for (int i = 0; i < *seed_size; i++) {
		points[seeds[i]].cluster_id = cluster_id;
	}

	/* Remove the core point from seeds */
	for (int i = 0; i < *seed_size; i++) {
		if (seeds[i] == point_idx) {
			seeds[i] = seeds[--(*seed_size)];
			break;
		}
	}

	/* Process all seed points */
	int current_seed = 0;
	while (current_seed < *seed_size) {
		int current_point = seeds[current_seed];

		/* Find neighbors of current seed point */
		int neighbor_count = cdbscan_region_query_custom(
			points, num_points, current_point, params, neighbors);

		if (neighbor_count >= params->min_pts) {
			/* Current point is also a core point */
			for (int i = 0; i < neighbor_count; i++) {
				int neighbor_idx = neighbors[i];

				if (points[neighbor_idx].cluster_id ==
					    CDBSCAN_UNCLASSIFIED ||
				    points[neighbor_idx].cluster_id ==
					    CDBSCAN_NOISE) {
					if (points[neighbor_idx].cluster_id ==
					    CDBSCAN_UNCLASSIFIED) {
						/* Add to seeds if it was unclassified */
						seeds[(*seed_size)++] =
							neighbor_idx;
					}

					/* Assign to current cluster */
					points[neighbor_idx].cluster_id =
						cluster_id;
				}
			}
		}

		current_seed++;
	}

	return 1; /* Successfully expanded cluster */
}

/* Expand cluster from a core point using KD-tree */
static int expand_cluster_kdtree(cdbscan_point_t *points, int num_points,
				 int point_idx, int cluster_id,
				 const cdbscan_params_t *params,
				 const kdtree_t *tree, int *neighbors,
				 int *seeds, int *seed_size)
{
	/* Get initial seeds from KD-tree range query */
	*seed_size = kdtree_range_query(tree, point_idx, params->eps, seeds);

	if (*seed_size < params->min_pts) {
		/* Not a core point */
		points[point_idx].cluster_id = CDBSCAN_NOISE;
		return 0;
	}

	/* Assign cluster ID to all points in seeds */
	for (int i = 0; i < *seed_size; i++) {
		points[seeds[i]].cluster_id = cluster_id;
	}

	/* Remove the core point from seeds */
	for (int i = 0; i < *seed_size; i++) {
		if (seeds[i] == point_idx) {
			seeds[i] = seeds[--(*seed_size)];
			break;
		}
	}

	/* Process all seed points */
	int current_seed = 0;
	while (current_seed < *seed_size) {
		int current_point = seeds[current_seed];

		/* Find neighbors of current seed point using KD-tree */
		int neighbor_count = kdtree_range_query(tree, current_point,
							params->eps, neighbors);

		if (neighbor_count >= params->min_pts) {
			/* Current point is also a core point */
			for (int i = 0; i < neighbor_count; i++) {
				int neighbor_idx = neighbors[i];

				if (points[neighbor_idx].cluster_id ==
					    CDBSCAN_UNCLASSIFIED ||
				    points[neighbor_idx].cluster_id ==
					    CDBSCAN_NOISE) {
					if (points[neighbor_idx].cluster_id ==
					    CDBSCAN_UNCLASSIFIED) {
						/* Add to seeds if it was unclassified */
						seeds[(*seed_size)++] =
							neighbor_idx;
					}

					/* Assign to current cluster */
					points[neighbor_idx].cluster_id =
						cluster_id;
				}
			}
		}

		current_seed++;
	}

	return 1; /* Successfully expanded cluster */
}

/* Utility functions */
cdbscan_point_t *cdbscan_create_points(int num_points, int dimensions)
{
	if (num_points <= 0 || dimensions <= 0) {
		return NULL;
	}

	cdbscan_point_t *points =
		(cdbscan_point_t *)calloc(num_points, sizeof(cdbscan_point_t));
	if (!points) {
		return NULL;
	}

	for (int i = 0; i < num_points; i++) {
		points[i].coords = (double *)calloc(dimensions, sizeof(double));
		if (!points[i].coords) {
			/* Clean up on failure */
			for (int j = 0; j < i; j++) {
				free(points[j].coords);
			}
			free(points);
			return NULL;
		}
		points[i].dimensions = dimensions;
		points[i].cluster_id = CDBSCAN_UNCLASSIFIED;
		points[i].index = i;
	}

	return points;
}

void cdbscan_free_points(cdbscan_point_t *points)
{
	/* This function assumes individual allocation per point */
	/* Users should free based on how they allocated */
	if (points) {
		free(points);
	}
}

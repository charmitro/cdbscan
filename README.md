# cdbscan

A C library implementing the DBSCAN (Density-Based Spatial Clustering of Applications with Noise) algorithm, based on the original KDD-96 paper by Ester et al.

## Building

```bash
$ make
$ make test                 # Run tests
$ make install              # Install to /usr/local
$ make install PREFIX=/path # Custom installation path
```

## Usage

```c
#include <cdbscan.h>

// Create points
cdbscan_point_t* points = cdbscan_create_points(num_points, 2);

// Populate with your data
for (int i = 0; i < num_points; i++) {
	points[i].coords[0] = x_data[i];
	points[i].coords[1] = y_data[i];
}

// Set parameters
cdbscan_params_t params = {
	.eps = 0.5,                         // Neighborhood radius
	.min_pts = 4,                       // Minimum points for core
	.dist_type = CDBSCAN_DIST_EUCLIDEAN,
	.use_kdtree = 1                     // Enable O(n log n) mode
};

// Run clustering
int num_clusters = cdbscan_cluster(points, num_points, params);

// Check results (cluster_id: >=0 cluster, -2 noise)
for (int i = 0; i < num_points; i++) {
	printf("Point %d: cluster %d\n", i, points[i].cluster_id);
}
```

## Examples

```bash
$ make examples
$ ./examples/example                 # Basic usage
$ ./examples/example_kdtree          # Performance comparison
```

## Authors

* Charalampos Mitrodimas <charmitro@posteo.net>
* The cdbscan developers

## License

GNU General Public License v3.0 - see LICENSE file for details.

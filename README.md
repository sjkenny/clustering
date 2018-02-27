# clustering
### DBSCAN and flood-fill clustering

Density- and intensity- based clustering methods enable feature extraction of point clouds and images, respectively, without *a priori* knowledge of the number of clusters.  These methods also allow for the identification of clusters with arbitrary sizes and shapes, as long as the density across the area of the cluster is consistent.  They are also tolerant to noise and outliers. As such, they are useful tools for feature extraction and quantification in imaging.

### Dependencies
[STORM image I/O](https://github.com/sjkenny/common) for STORM molecule lists generated with Insight3
- OpenMolListTxt for .txt molecule lists
- OpenMolList for .bin molecule lists
- WriteMolBinNXcYcZc for saving molecule list indexed by cluster

### Usage

- Call dbscan_fcn.m and specify x,y point cloud coordinates as column vectors, as well as MinPts and Eps for minimum density and search radius, respectively.  
- Call flood_fill_indices.m with a grayscale or binary image input.  For grayscale images, nonzero pixels will be connected.
- Both functions return indexed lists corresponding to the molecules or pixels in each cluster, as well as centers-of-mass and sizes (either # of points or area) of each cluster.

#### DBSCAN

[Density-based spatial clustering of applications with noise (DBSCAN)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.121.9220) is a point cloud clustering algorithm which clusters points with sufficent spatial density.  This implementation takes two parameters: a search radius (Eps) and a minimum number of points (MinPts).  Clustering proceeds by first finding all neighbors within Eps of a point; if the number of neighbors is above the threshold MinPts, these points are added to a stack, and then queried in the same way.  If a point has fewer than MinPts neighbors, it is considered a cluster edge, and its neighbors are not queried.  Once all points in the stack have been queried, the cluster is complete, and the next point in the point cloud initializes a new potential cluster.

Shown here is the (accelerated) search progression of DBSCAN for a small point cloud with two clusters.  Points with sufficient local density are colored magenta, and are added to the cluster.  Once the first cluster is finished, the second is initialized automatically.

![DBSCAN clustering](https://i.imgur.com/rRei4Mq.gif)

DBSCAN accurately identifies clusters of disparate sizes and shapes:

![DBSCAN clustering](https://i.imgur.com/Q6aqOb9.gif)

#### Flood-fill

Flood-fill is an algorithm which connects neighboring pixels of the same value.  It can be applied directly to grayscale image data, or to point cloud data that has been binned into pixels.  A minimum intensity threshold for connection must be specified, or the input image must be cast to binary.  This implementation automatically connects all regions of like-valued pixels, and returns an indexed output matrix of pixels belonging to each cluster.  It can be used as a fast alternative to DBSCAN for very large point clouds, when some pixelation noise is tolerable.

This implementation is a queue-based search along the Y direction, first finding the minimum and maximum X values of connected pixels, and then filling in the result row-by-row, as shown here:

![Flood-fill clustering](https://i.imgur.com/9WJZccO.gif)

# IIT-Patna-intern
IIT-Patna internship on the topic k-means clustering algorithm...implemented in C
The Algorithm is explained as follows:

k-means clustering, or Lloyd's algorithm [2], is an iterative, data-partitioning algorithm that assigns n observations to exactly one of k clusters defined by centroids, where k is chosen before the algorithm starts.

The algorithm proceeds as follows:

1. Choose k initial cluster centers (centroid). For example, choose k observations at random (by using 'Start','sample') or use the k-means ++ algorithm for cluster center initialization (the default).

2. Compute point-to-cluster-centroid distances of all observations to each centroid.

3. There are two ways to proceed (specified by OnlinePhase):

	1. Batch update — Assign each observation to the cluster with the closest centroid.

	2. Online update — Individually assign observations to a different centroid if the reassignment decreases the sum of the within-cluster, sum-of-squares point-to-cluster-centroid distances.

	For more details, see Algorithms.

4. Compute the average of the observations in each cluster to obtain k new centroid locations.

5. Repeat steps 2 through 4 until cluster assignments do not change, or the maximum number of iterations is reached.

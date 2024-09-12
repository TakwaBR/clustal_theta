"""
This module performs hierarchical clustering using the UPGMA.

Dependencies
------------
- pandas: For handling the distance matrix and data manipulation.
- numpy: For numerical operations, especially for handling NaN
    values and distance computations.
- matplotlib: For plotting the dendrogram.
- scipy: For creating the dendrogram.

Functions
---------
- compute_distance(dist1: float, dist2: float, size1: int, size2: int) -> float:
    Calculate the weighted average distance between two clusters based on their
    distances to a third cluster and their sizes.
    
- hierarchical_clustering(dist_matrix: pd.DataFrame) -> Tuple[List[List[float]],
str]:
    Perform hierarchical clustering using the UPGMA method on a distance
    matrix. Returns a linkage matrix and the final cluster formed.

Example
-------
To use this module, you need to provide a distance matrix as a pandas
DataFrame. The script includes an example distance matrix and generates
a dendrogram plot showing the hierarchical clustering result.

Usage
-----
1. Define your distance matrix as a dictionary and convert it to a
pandas DataFrame.
2. Call the `hierarchical_clustering` function with this DataFrame
to perform clustering.
3. The resulting linkage matrix is printed, and a dendrogram is plotted 
showing the clustering result.

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram


def compute_distance(dist1, dist2, size1, size2):
    """
    Calculate the weighted average distance between two clusters.

    Parameters
    ----------
    dist1 : float
        The distance from the first cluster to the third cluster.
    dist2 : float
        The distance from the second cluster to the third cluster.
    size1 : int
        The size (number of elements) of the first cluster.
    size2 : int
        The size (number of elements) of the second cluster.

    Returns
    -------
    float
        The computed distance between the two clusters.
    """
    return (size1 * dist1 + size2 * dist2) / (size1 + size2)


def hierarchical_clustering(dist_matrix):
    """
    Perform hierarchical clustering using the UPGMA method.

    Parameters
    ----------
    dist_matrix : pd.DataFrame
        A DataFrame representing the distance matrix where the rows
        and columns are sequence names and the cell values are distances.

    Returns
    -------
    Tuple[List[List[float]], str]
        A tuple where the first element is a list of lists representing
        the linkage matrix and the second element is the final cluster formed.
    """
    cluster_sizes = {seq: 1 for seq in dist_matrix.columns}
    linkage_matrix = []
    cluster_indices = {seq: idx for idx, seq in enumerate(dist_matrix.columns)}

    while len(dist_matrix) > 1:
        # Searching for the min distance
        min_dist = dist_matrix.min().min()
        min_loc = np.where(dist_matrix == min_dist)
        row_idx, col_idx = min_loc[0][0], min_loc[1][0]
        row_label = dist_matrix.index[row_idx]
        col_label = dist_matrix.columns[col_idx]

        # Creation of a new cluster and a new row for the cluster
        new_cluster = (row_label, col_label)
        new_row = {new_cluster: [np.nan]}
        remaining_seqs = [seq for seq in dist_matrix.columns if seq not in
                          new_cluster]

    for seq in remaining_seqs:
        # Check if the sequence is to the right of both clusters in the matrix
        if (row_idx < dist_matrix.columns.get_loc(seq)
            and col_idx < dist_matrix.columns.get_loc(seq)):
            # Calculate the distance from the new cluster to the sequence
            # based on distances to both clusters
            distance = compute_distance(dist_matrix.iloc[row_idx, dist_matrix.columns.get_loc(seq)], 
                                        dist_matrix.iloc[col_idx, dist_matrix.columns.get_loc(seq)], 
                                        cluster_sizes[row_label], cluster_sizes[col_label])
        # Check if the sequence is to the right of the first cluster but
        # to the left of the second cluster
        elif row_idx < dist_matrix.columns.get_loc(seq) and col_idx > dist_matrix.columns.get_loc(seq):
            # Calculate the distance from the new cluster to the
            # sequence based on distance to the first cluster (row) and second
            # cluster (column)
            distance = compute_distance(dist_matrix.iloc[row_idx, dist_matrix.columns.get_loc(seq)], 
                                        dist_matrix.iloc[dist_matrix.columns.get_loc(seq), col_idx], 
                                        cluster_sizes[row_label], cluster_sizes[col_label])
        # Check if the sequence is to the left of the first cluster
        # but to the right of the second cluster
        elif row_idx > dist_matrix.columns.get_loc(seq) and col_idx < dist_matrix.columns.get_loc(seq):
            # Calculate the distance from the new cluster to the
            # sequence based on distance to the second cluster (row) and
            # first cluster (column)
            distance = compute_distance(dist_matrix.iloc[dist_matrix.columns.get_loc(seq), row_idx], 
                                        dist_matrix.iloc[col_idx, dist_matrix.columns.get_loc(seq)], 
                                        cluster_sizes[row_label], cluster_sizes[col_label])
        else:
            # The sequence is to the left of both clusters in the matrix
            # Calculate the distance from the new cluster to the sequence
            # based on distances to both clusters
            distance = compute_distance(dist_matrix.iloc[dist_matrix.columns.get_loc(seq), row_idx], 
                                        dist_matrix.iloc[dist_matrix.columns.get_loc(seq), col_idx], 
                                        cluster_sizes[row_label], cluster_sizes[col_label])

        # Append the computed distance to the new row for the new cluster
        new_row[new_cluster].append(distance)


        # Update linkage matrix
        linkage_matrix.append([cluster_indices[row_label],
                               cluster_indices[col_label], min_dist,
                               cluster_sizes[row_label] +
                               cluster_sizes[col_label]])

        # Transform the matrix
        # Drop the two grouped sequences from rows and columns
        dist_matrix.drop([row_label, col_label], axis=0, inplace=True)
        dist_matrix.drop([row_label, col_label], axis=1, inplace=True)

        # Update cluster sizes
        cluster_sizes[new_cluster] = (cluster_sizes.pop(row_label) +
                                      cluster_sizes.pop(col_label))
        cluster_indices[new_cluster] = len(cluster_indices)

        # Add the new row and column
        dist_matrix.insert(0, new_cluster, [np.nan]*len(dist_matrix))

        new_row_df = pd.DataFrame([new_row[new_cluster]],
                                  columns=dist_matrix.columns,
                                  index=[new_cluster])
        dist_matrix = pd.concat([new_row_df, dist_matrix])

    return linkage_matrix, new_cluster


if __name__ == "__main__":
    DIST_MATRIX = {'Hu': [np.nan, np.nan, np.nan, np.nan, np.nan],
                   'Ch': [15, np.nan, np.nan, np.nan, np.nan],
                   'Go': [45, 30, np.nan, np.nan, np.nan],
                   'Or': [143, 126, 92, np.nan, np.nan],
                   'Gi': [198, 179, 179, 179, np.nan]}

    test_df = pd.DataFrame(DIST_MATRIX, index=['Hu', 'Ch', 'Go', 'Or', 'Gi'])
    linkage_result, clusters = hierarchical_clustering(test_df)
    print(linkage_result)
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_result, labels=['Hu', 'Ch', 'Go', 'Or', 'Gi'])
    plt.title("Dendrogram")
    plt.xlabel("Clusters")
    plt.ylabel("Distance")
    plt.show()

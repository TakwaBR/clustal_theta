import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

def compute_distance(dist1, dist2, size1, size2):
    return (size1 * dist1 + size2 * dist2) / (size1 + size2)

def hierarchical_clustering(dist_matrix):
    cluster_sizes = {seq: 1 for seq in dist_matrix.columns}
    linkage_matrix = []
    cluster_indices = {seq: idx for idx, seq in enumerate(dist_matrix.columns)}
    
    while len(dist_matrix) > 1:
        min_dist = dist_matrix.min().min()
        min_loc = np.where(dist_matrix == min_dist)
        row_idx, col_idx = min_loc[0][0], min_loc[1][0]
        row_label = dist_matrix.index[row_idx]
        col_label = dist_matrix.columns[col_idx]

        new_cluster = (row_label, col_label)
        new_row = {new_cluster: [np.nan]}
        remaining_seqs = [seq for seq in dist_matrix.columns if seq not in new_cluster]
        
        for seq in remaining_seqs:
            if row_idx < dist_matrix.columns.get_loc(seq) and col_idx < dist_matrix.columns.get_loc(seq):
                distance = compute_distance(dist_matrix.iloc[row_idx, dist_matrix.columns.get_loc(seq)], 
                                            dist_matrix.iloc[col_idx, dist_matrix.columns.get_loc(seq)], 
                                            cluster_sizes[row_label], cluster_sizes[col_label])
            elif row_idx < dist_matrix.columns.get_loc(seq) and col_idx > dist_matrix.columns.get_loc(seq):
                distance = compute_distance(dist_matrix.iloc[row_idx, dist_matrix.columns.get_loc(seq)], 
                                            dist_matrix.iloc[dist_matrix.columns.get_loc(seq), col_idx], 
                                            cluster_sizes[row_label], cluster_sizes[col_label])
            elif row_idx > dist_matrix.columns.get_loc(seq) and col_idx < dist_matrix.columns.get_loc(seq):
                distance = compute_distance(dist_matrix.iloc[dist_matrix.columns.get_loc(seq), row_idx], 
                                            dist_matrix.iloc[col_idx, dist_matrix.columns.get_loc(seq)], 
                                            cluster_sizes[row_label], cluster_sizes[col_label])
            else:
                distance = compute_distance(dist_matrix.iloc[dist_matrix.columns.get_loc(seq), row_idx], 
                                            dist_matrix.iloc[dist_matrix.columns.get_loc(seq), col_idx], 
                                            cluster_sizes[row_label], cluster_sizes[col_label])

            new_row[new_cluster].append(distance)
        
        # Update linkage matrix
        linkage_matrix.append([cluster_indices[row_label], cluster_indices[col_label], min_dist, cluster_sizes[row_label] + cluster_sizes[col_label]])
        
        # Transform the matrix
        # Drop the two grouped sequences from rows and columns
        dist_matrix.drop([row_label, col_label], axis=0, inplace=True)
        dist_matrix.drop([row_label, col_label], axis=1, inplace=True)
        
        # Update cluster sizes
        cluster_sizes[new_cluster] = cluster_sizes.pop(row_label) + cluster_sizes.pop(col_label)
        cluster_indices[new_cluster] = len(cluster_indices)
        
        # Add the new row and column
        dist_matrix.insert(0, new_cluster, [np.nan]*len(dist_matrix))
        
        new_row_df = pd.DataFrame([new_row[new_cluster]], columns=dist_matrix.columns, index=[new_cluster])
        dist_matrix = pd.concat([new_row_df, dist_matrix])

    return linkage_matrix

if __name__ == "__main__":
    DIST_MATRIX = {'Hu': [np.nan, np.nan, np.nan, np.nan, np.nan],
                   'Ch': [15, np.nan, np.nan, np.nan, np.nan],
                   'Go': [45, 30, np.nan, np.nan, np.nan],
                   'Or': [143, 126, 92, np.nan, np.nan],
                   'Gi': [198, 179, 179, 179, np.nan]}
    
    test_df = pd.DataFrame(DIST_MATRIX, index=['Hu', 'Ch', 'Go', 'Or', 'Gi'])
    linkage_result = hierarchical_clustering(test_df)
    
    # Plot the dendrogram
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_result, labels=['Hu', 'Ch', 'Go', 'Or', 'Gi'])
    plt.title("Dendrogram")
    plt.xlabel("Clusters")
    plt.ylabel("Distance")
    plt.show()

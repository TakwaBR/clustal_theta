"""
Module for sequence clustering and alignment.

This module provides functions to extract sequences from clusters
and perform pairwise alignments within specified clusters.
The primary functions included are:

- `extract_sequences_from_clusters`: Extracts sequences for each cluster
based on cluster IDs from a given collection of sequences.
- `align_seq`: Aligns sequences within a specified cluster using pairwise
alignment.

Note
----
This module is a work in progress. The functionality may be incomplete
or subject to change. Specifically, the alignment function is currently
set up to handle only 'Cluster1' and may require additional implementation
to handle other clusters or perform further operations.

Dependencies
------------
- `numpy`
- `pandas`
- `pair_alignment` (custom module for pairwise alignment)
- `upgma` (custom module for hierarchical clustering)

Example Usage
-------------
```python
cluster_structure = {
    'Cluster1': ['Homo', 'Gorilla'],
    'Cluster2': ['Nomascus', 'Hylobates'],
    'Cluster3': ['Trachypithecus', 'Chlorocebus'],
    # Add more clusters as needed
}

SEQUENCES = pa.SEQUENCES  # Custom sequence records loaded from pair_alignment
module

cluster_sequences = extract_sequences_from_clusters(clusters, SEQUENCES)
align_seq(cluster_sequences)
"""
import numpy as np
import pandas as pd
import pair_alignment as pa
import upgma as up


def extract_sequences_from_clusters(cluster_structure, sequences):
    """
    Extract sequences for each cluster from a collection of sequences.

    Parameters
    ----------
    cluster_structure : dict
        A dictionary where keys are cluster names and values are lists of
        sequence IDs.
    sequences : list
        A list of sequence records, where each record has an `id` and
        a `seq` attribute.

    Returns
    -------
    dict
        A dictionary where keys are cluster names and values are
        dictionaries of sequence IDs to sequences.
    """
    cluster_sequences = {}

    for cluster, ids in cluster_structure.items():
        cluster_sequences[cluster] = {
            seq_id: seq_record.seq
            for seq_record in sequences
            for seq_id in ids if seq_record.id == seq_id
        }

    return cluster_sequences


def align_seq(clusters_sequences):
    """
    Align sequences within a specified cluster using pairwise alignment.

    Parameters
    ----------
    clusters_sequences : dict
        A dictionary where keys are cluster names and values are dictionaries
        of sequence IDs to sequences.

    Returns
    -------
    None
        The function prints the alignment results for 'Cluster1'
        if it exists and contains at least two sequences.
    """
    aligned_seq = {}
    if 'Cluster1' in clusters_sequences:
        sequences = clusters_sequences['Cluster1']

        if len(sequences) >= 2:
            seq_ids = list(sequences.keys())[:2]
            seq1_index = next(i for i, seq in enumerate(SEQUENCES) if
                              seq.id == seq_ids[0])
            seq2_index = next(i for i, seq in enumerate(SEQUENCES) if
                              seq.id == seq_ids[1])

            print(f"Alignment for the cluster Cluster1:\n")
            pa.pair_alignment(seq1_index, seq2_index)


if __name__ == "__main__":
    cluster_structure = ((((('Ho', 'Go'), ('No', 'Hy')), (
        (('Tr', 'Ch'), 'Pi'), 'Ma')), 'Po'), 'Rh')

    clusters = {
        'Cluster1': ['Homo', 'Gorilla'],
        'Cluster2': ['Nomascus', 'Hylobates'],
        'Cluster3': ['Trachypithecus', 'Chlorocebus'],
        'Cluster4': ['Piliocolobus'],
        'Cluster5': ['Macaca'],
        'Cluster6': ['Pongo'],
        'Cluster7': ['Rhinopithecus']
    }
    SEQUENCES = pa.SEQUENCES
    cluster_sequences = extract_sequences_from_clusters(clusters, SEQUENCES)
    align_seq(cluster_sequences)

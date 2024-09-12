"""
Sequence Alignment and Distance Calculation Module.

This module performs pairwise sequence alignment for a set of sequences,
calculates the resulting scores, and computes a distance matrix based on
these scores. The distance matrix is then used for further clustering 
analysis using the UPGMA algorithm.

Dependencies:
- numpy
- pandas
- Bio (Biopython)
- re (regular expressions)
- cProfile (for profiling)
- upgma (for hierarchical clustering)

Global Variables:
- FASTA_FILE (str): Path to the FASTA file containing sequences.
- BLOSUM_FILE (str): Path to the CSV file containing the BLOSUM62 matrix.
- SEQUENCES (list): List of sequences parsed from the FASTA file.
- BLOSUM_MATRIX (pandas.DataFrame): DataFrame containing the BLOSUM62 matrix.
- GAP_PENALTY (int): Gap penalty used in sequence alignment.
- AA_INDEX (dict): Dictionary mapping amino acid characters to their indices in the BLOSUM62 matrix.

Functions:
- matrix_init(sequence1_index, sequence2_index):
    Initializes the score matrix for sequence alignment with gap penalties.
    
- display_alignment(aligned_seq1, aligned_seq2, block_size=80):
    Displays the aligned sequences with alignment bars in blocks.
    
- pair_alignment(sequence1_index, sequence2_index):
    Performs pairwise sequence alignment without an explicit traceback matrix.
    
- calculate_distance(all_score_matrix):
    Calculates the distance matrix from a given score matrix.
    
- align_all_pairs():
    Performs pairwise sequence alignment for all pairs of sequences, calculates the alignment scores, and computes a distance matrix.

Usage:
- The module is designed to be run as a script. It aligns all sequences, calculates the distance matrix, and applies the UPGMA clustering algorithm.
- Profiling functions for performance analysis are provided but commented out by default.

Example:
    To run the module and perform clustering, execute the script. It will align all sequences, calculate distances, and perform UPGMA clustering.

    >>> python this_script.py
"""
import numpy as np
import pandas as pd
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import upgma

# Loading the fasta file and the Blosum62 matrix
FASTA_FILE = "../data/p53_sequences.fasta"
BLOSUM_FILE = "../data/blosum62.csv"
# Reading the sequences in the fasta file
SEQUENCES = list(SeqIO.parse(FASTA_FILE, "fasta"))
BLOSUM_MATRIX = pd.read_csv(BLOSUM_FILE, index_col=0)
# Choosing the gap penalty
GAP_PENALTY = -10

np.set_printoptions(precision=3, suppress=True)

# Pre-calculation of indices for amino acids
AA_INDEX = {aa: i for i, aa in enumerate(BLOSUM_MATRIX.index)}


def matrix_init(sequence1_index, sequence2_index):
    """
    Initialize the score matrix for sequence alignment with gap penalties.

    Args:
        sequence1_index (int): Index of the first sequence in the SEQUENCES
        list.
        sequence2_index (int): Index of the second sequence in the SEQUENCES
        list.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: The initialized score matrix with gap penalties.
            - str: The first sequence corresponding to the given index.
            - str: The second sequence corresponding to the given index.

    Notes:
        - The score matrix is initialized with dimensions (len(sequence1) + 1,
        len(sequence2) + 1).
        - The first row and column of the matrix are initialized with
        cumulative gap penalties.
        - The sequences are retrieved from a global SEQUENCES list using the
        provided indices.
    """
    sequence1 = str(SEQUENCES[sequence1_index].seq)
    sequence2 = str(SEQUENCES[sequence2_index].seq)

    score_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1))

    # Initialization of the first row and column with gap penalties
    score_matrix[0, :] = np.arange(0, (len(sequence2) + 1) * GAP_PENALTY,
                                   GAP_PENALTY)
    score_matrix[:, 0] = np.arange(0, (len(sequence1) + 1) * GAP_PENALTY,
                                   GAP_PENALTY)
    return score_matrix, sequence1, sequence2


def display_alignment(aligned_seq1, aligned_seq2, block_size=80):
    """
    Display the aligned sequences with alignment bars in blocks.

    Args:
        aligned_seq1 (str): The first aligned sequence.
        aligned_seq2 (str): The second aligned sequence.
        block_size (int, optional): The number of characters per block for
        display. Defaults to 80.

    Returns:
        None

    Notes:
        - The function prints the aligned sequences in blocks of `block_size`
        characters.
        - Alignment bars (`|` for matches and `*` for mismatches) are displayed
        between the sequences.
        - The alignment bars indicate the positional alignment of characters
        between the two sequences.

    Example:
        >>> display_alignment("AGCTAGCTAG", "AGCTTGCTAG")
        AGCTAGCTAG
        |||*|*||||
        AGCTTGCTAG
    """
    for start in range(0, len(aligned_seq1), block_size):
        block_seq1 = aligned_seq1[start:start + block_size]
        block_seq2 = aligned_seq2[start:start + block_size]

        align_line = []
        for a1, a2 in zip(block_seq1, block_seq2):
            if a1 == a2:
                align_line.append('|')
            else:
                align_line.append('*')

        print("".join(block_seq1))
        print("".join(align_line))
        print("".join(block_seq2))
        print()


def pair_alignment(sequence1_index, sequence2_index):
    """
    Perform pairwise sequence alignment without an explicit traceback matrix.

    Args:
        sequence1_index (int): Index of the first sequence in the SEQUENCES
        list.
        sequence2_index (int): Index of the second sequence in the SEQUENCES
        list.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: The score matrix after alignment calculations.
            - list: The aligned first sequence.
            - list: The aligned second sequence.

    Notes:
        - The function initializes a score matrix with gap penalties using
        the `matrix_init` function.
        - The sequences are aligned using scores from the BLOSUM matrix and
        gap penalties.
        - The alignment is obtained by backtracking through the score matrix.
        - The resulting aligned sequences are returned in their original order.
        - The sequences are represented as lists of characters for the
        alignment output.
    """
    # Matrix initialization
    score_matrix, seq1, seq2 = matrix_init(sequence1_index, sequence2_index)

    # Convert amino acid sequences to indices for accessing BLOSUM scores
    seq1_indices = []
    for aa in seq1:
        seq1_indices.append(AA_INDEX[aa])

    seq2_indices = []
    for aa in seq2:
        seq2_indices.append(AA_INDEX[aa])

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score = BLOSUM_MATRIX.iloc[seq1_indices[i - 1],
                                       seq2_indices[j - 1]]

            match = score_matrix[i - 1, j - 1] + score
            delete = score_matrix[i - 1, j] + GAP_PENALTY
            insert = score_matrix[i, j - 1] + GAP_PENALTY
            # Choose the maximum score among match, delete, and insert
            score_matrix[i, j] = max(match, delete, insert)

    # Backtraking
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        current_score = score_matrix[i, j]
        # Check if the current score is from a match
        if (i > 0 and j > 0 and current_score == score_matrix[i - 1, j - 1]
                + BLOSUM_MATRIX.iloc[seq1_indices[i - 1],
                                     seq2_indices[j - 1]]):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        # Check if the current score is from a deletion
        elif i > 0 and current_score == score_matrix[i - 1, j] + GAP_PENALTY:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        # The current score is from an insertion
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
    # Reverse the aligned sequences to get the correct order
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    display_alignment(aligned_seq1, aligned_seq2)

    return score_matrix, aligned_seq1, aligned_seq2


def calculate_distance(all_score_matrix):
    """
    Calculate the distance matrix from a given score matrix.

    Args:
        all_score_matrix (numpy.ndarray): The score matrix used to calculate
        distances.

    Returns:
        numpy.ndarray: The distance matrix, where distances are scaled between
        0 and 1.
    """
    distance_matrix = np.zeros((len(all_score_matrix), len(all_score_matrix)))
    max_score = np.max(all_score_matrix)
    min_score = np.min(all_score_matrix)
    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix)):
            distance_matrix[i, j] = (1 - (all_score_matrix[i, j]-min_score) /
                                     (max_score-min_score))

    return distance_matrix


def align_all_pairs():
    """
    Perform pairwise sequence alignment for all pairs of sequences.

    This function aligns all pairs of sequences using the pairwise alignment
    method, calculates the alignment scores,
    and then computes a distance matrix based on these scores.
    The distance matrix is scaled to a [0, 1] range, with higher
    scores representing closer similarity and lower distances.
    The distance matrix is then returned as a DataFrame.

    Returns:
        pandas.DataFrame: A DataFrame representing the distance matrix, where
        rows and columns are indexed by the organism abbreviations.
    """
    all_score_matrix = np.zeros((len(SEQUENCES), len(SEQUENCES)))
    organismes = []
    for i in range(len(SEQUENCES)):
        # Extract organism abbreviation from the sequence description
        description_i = SEQUENCES[i].description
        organism_match_i = re.search(r"\[(.*?)\]", description_i)
        organism_i = organism_match_i.group(1)
        organismes.append(organism_i[0:2])
        for j in range(i + 1, len(SEQUENCES)):
            # Perform pairwise alignment for sequences i and j
            print(f"#######################################################\n")
            description_j = SEQUENCES[j].description
            organism_match_j = re.search(r"\[(.*?)\]", description_j)
            organism_j = organism_match_j.group(1)
            print(f"Pair alignment between {organism_i} and {organism_j}\n")
            score_matrix, aligned_seq1, aligned_seq2 = pair_alignment(i, j)
            # Store the final score from the score matrix in the score_matrix
            all_score_matrix[i, j] = score_matrix[-1, -1]

    distance_matrix = calculate_distance(all_score_matrix)
    distance_matrix[distance_matrix == 1] = np.nan

    # Create a DataFrame for better visualization and usability
    distance_df = pd.DataFrame(
        distance_matrix,
        columns=organismes,
        index=organismes
    )
    return distance_df


if __name__ == "__main__":
    distance_matrix = align_all_pairs()
    linkage_mat = upgma.hierarchical_clustering(distance_matrix)

    plt.figure(figsize=(10, 7))
    dendrogram(linkage_mat, labels=['Ho', 'Po', 'Go', 'No', 'Hy', 'Rh', 'Tr', 'Ch', 'Pi', 'Ma'])
    plt.title("Dendrogramme")
    plt.xlabel("Clusters")
    plt.ylabel("Distance")
    plt.show()

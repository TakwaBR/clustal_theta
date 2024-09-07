import numpy as np
import pandas as pd
from Bio import SeqIO
import cProfile

FASTA_FILE = "../data/p53_sequences.fasta"
BLOSUM_FILE = "../data/blosum62.csv"
SEQUENCES = list(SeqIO.parse(FASTA_FILE, "fasta"))
BLOSUM_MATRIX = pd.read_csv(BLOSUM_FILE)

gap_penalty = -10

def matrix_initialisation(sequence1_index, sequence2_index):
      
    sequence1 = str(SEQUENCES[sequence1_index].seq)
    sequence2 = str(SEQUENCES[sequence2_index].seq)
    
    score_matrix = pd.DataFrame(
        0,
        index= [""] + list(sequence1),
        columns= [""] + list(sequence2)
    )

    decrement_values_row1 = np.arange(10, 10 * len(score_matrix.columns) + 10, 10)
    decrement_values_col1 = np.arange(10, 10 * score_matrix.shape[0] + 10, 10)
    score_matrix.iloc[0] -= decrement_values_row1
    score_matrix.iloc[1:,0] -= decrement_values_col1[1:]

    trace_matrix = pd.DataFrame("", index=list(sequence1), columns=list(sequence2))
    trace_matrix.iloc[0, 1:] = "L"  # Left for first row
    trace_matrix.iloc[1:, 0] = "U"  # Up for first column

    return score_matrix, sequence1, sequence2


def pair_alignement(sequence1_index, sequence2_index):
    score_matrix, seq1, seq2 = matrix_initialisation(sequence1_index, sequence2_index)
    
    # Calculer les scores
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # Obtenir les acides aminés
            aa1 = seq1[i - 1]
            aa2 = seq2[j - 1]

            score = BLOSUM_MATRIX.loc[aa1, aa2]

            # Calculer les valeurs possibles pour la cellule actuelle
            match = score_matrix.iloc[i - 1, j - 1] + score
            delete = score_matrix.iloc[i - 1, j] + gap_penalty
            insert = score_matrix.iloc[i, j - 1] + gap_penalty

            # Choisir la meilleure option
            score_matrix.iloc[i, j] = max(match, delete, insert)
    
    #print(score_matrix)


if __name__ == "__main__" :
    #pair_alignement(1, 2)
    #cProfile.run('pair_alignement(1, 2)')
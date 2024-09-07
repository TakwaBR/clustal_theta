import numpy as np
import pandas as pd
from Bio import SeqIO

FASTA_FILE = "../data/p53_sequences.fasta"
SCORES_FILE = "../data/blosum62.csv"
SEQUENCES = list(SeqIO.parse(FASTA_FILE, "fasta"))
SCORE_MATRIX = pd.read_csv(SCORES_FILE)
print(SCORE_MATRIX)
#def pair_alignement(sequence1, sequence2):

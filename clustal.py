import numpy as np
import pandas as pd
from Bio import SeqIO

FASTA_FILE = "../data/sequences.fasta"

sequences = list(SeqIO.parse(FASTA_FILE, "fasta"))



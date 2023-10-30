import pandas as pd
from seq_mod import *
from Bio import Align

# read the fasta file first
fasta_name = "mouse.fasta"
with open(fasta_name) as fin:
    # Split the protein entries with ">" in the beginning of the description
    split_list = [protein for protein in fin.read()[1:].split('\n>')]
    # get the gene name by splitting the data with "="
    gene_list = [entry.split("=")[3].split(" ")[0] for entry in split_list]
    # Protein sequence with join for multiple lines
    protein_list = ["".join(split_list[uniport].split("\n")[1:]) for uniport in range(0, len(split_list))]

# Read the xlsx file of mutaed segments, only keep the useful columns
# Also drop the rows where there is no mutation sequence
mutated = (pd.read_excel("Neoantigen and Damps_zz.xlsx", usecols=[1, 2, 3]).dropna(subset="Neopeptide sequence")
           .fillna("No_data"))

for index, row in mutated.iterrows():
    if row["Gene origin"] != "No_data":
        row["ID"] = row["Gene origin"]
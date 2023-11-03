import pandas as pd
from seq_mod import *
import numpy as np
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
mutated = (pd.read_excel("Neoantigen and Damps_zz.xlsx", usecols=[0, 1, 2, 3]).dropna(subset="Neopeptide sequence")
           .fillna("No_data"))
# Modify the data
for index, row in mutated.iterrows():
    if row["Gene origin"] != "No_data":
        row["ID"] = row["Gene origin"]
mutated = mutated[mutated["Cell line"] != "DAMP"].iloc[:, 1:3]

# Long string for aho-corasick
aho_seq = "|".join(protein_list)
zero_list = np.zeros(len(aho_seq), dtype= int).tolist()
ind = int(0)
for i in range(0,len(zero_list)):
    if aho_seq[i] != "|":
        zero_list[i] += ind
    else:
        ind += 1


for index, mutation in mutated.iterrows():
        if "[" in mutation["Neopeptide sequence"]:
            original = search_seq(mutation["Neopeptide sequence"])[0]
            modified = search_seq(mutation["Neopeptide sequence"])[1]
            # Locate the sequence with gene name, some don't have gene name matched
            try:
                seq_index = gene_list.index(mutation["ID"])
                seq_string = protein_list[seq_index]
                # replace the original string with the modified string
                final_sequence = seq_string.replace(original, modified)
                name = split_list[seq_index].split("\n")[0]
                # Get the result in fasta format
                fasta_protein_write = (
                    (">" + name.split("|")[0] + "|PPP{}|" + name.split("|")[2] + "\n" + final_sequence)
                    .format(str(str(index).zfill(3)) + "_" + str(1)))
            except ValueError:
                n = 0
                for protein in protein_list:
                    if protein.find(original) != -1:
                        n += 1
                        final_sequence = protein.replace(original, modified)
                        name = split_list[protein_list.index(protein)].split("\n")[0]
                        fasta_protein_write = (
                            (">" + name.split("|")[0] + "|PPP{}|" + name.split("|")[2] + "\n" + final_sequence)
                            .format(str(str(index).zfill(3)) + "_" + str(n)))

        # The segment do not have the bracket, but have gene origin name
        else:
            continue
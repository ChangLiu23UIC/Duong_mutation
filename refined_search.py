
# This is meant to
def fine_search(seq_dict: dict, protein_list: list):
    for seq, protein in seq_dict:
        seq_list = [seq[i:i+5] for i in range(len(seq) - 5 + 1)]
        for i in seq_list:

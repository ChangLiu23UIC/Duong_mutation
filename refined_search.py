# This is meant to
def fine_search(seq_dict: dict, protein_list: list):
    for seq, protein in seq_dict.items():
        best = 0
        res_protein = 0
        seq_list = [seq[i:i + 5] for i in range(len(seq) - 5 + 1)]
        for sub_protein in protein:
            count = 0
            for i in seq_list:
                if i in protein_list[sub_protein]:
                    count += 1
            if count > best:
                best = count
                res_protein = sub_protein
        seq_dict[seq] = res_protein

    res_dict = {k: v for k, v in seq_dict.items() if v != 0}
    return res_dict

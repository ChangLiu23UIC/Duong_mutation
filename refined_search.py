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

def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

def write_append(seq_dict: dict, write_list: list, protein_list: list, split_list: list, protein_index: int):
    for sequence, protein in seq_dict.items():
        protein_seq = protein_list[protein]
        position = protein_seq.find(sequence[0:6])
        name = split_list[protein].split("\n")[0]
        final_sequence_1 = protein_seq[:position] + sequence + protein_seq[(position + len(sequence)):]
        final_sequence = insert_newlines(final_sequence_1)
        fasta_protein_write = (
            (">" + name.split("|")[0] + "|PPP{}|" + name.split("|")[2] + "\n" + final_sequence + "\n")
            .format(str(str(protein_index).zfill(3))))
        protein_index += 1
        write_list.append(fasta_protein_write)

    return protein_index, write_list
import re


# This is the function that takes in the sequence of mutation and return a list of string
# for searching and replacing [original,modified]
def search_seq(unmodified: str) -> list:
        # regular expression for sequence with brackets
        parameter = re.compile(r"\[[A-Z]\/[A-Z]\]")
        # Search with the parameter
        position = parameter.search(unmodified)
        # Generate the original sequence and the mutation sequence with the bracket index
        original = unmodified[0:position.start()] + unmodified[position.start() + 1] + unmodified[position.end():]
        modified = unmodified[0:position.start()] + unmodified[position.start() + 3] + unmodified[position.end():]
        return [original, modified]


# Do the search of the aho-corasick
"""for end_index, (insert_order, original_value) in automaton.iter(aho_seq):

    start_index = end_index - len(original_value) + 1
    assert aho_seq[start_index:start_index + len(original_value)] == original_value
    if original_value in seq_dict:
        seq_dict[no_brac_no_gene.iloc(0)[insert_order]["Neopeptide sequence"]].append(zero_list[end_index])
    else:
        seq_dict[no_brac_no_gene.iloc(0)[insert_order]["Neopeptide sequence"]] = []
        seq_dict[no_brac_no_gene.iloc(0)[insert_order]["Neopeptide sequence"]].append(zero_list[end_index])
        """


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

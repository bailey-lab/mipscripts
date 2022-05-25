from collections import Counter
import csv
import os
import re
from string import punctuation, whitespace

def header_to_snake_case(path, overwrite=True):
    """Converts header column names of a file to snake case.

    Args:
        path (str): The path to the file.
        overwrite (bool): A boolean indicating whether to overwrite the original
            file. If `False`, a new file will be generated with the same name
            as the original file but with `_fix` appended to it.
    """
    with open(path) as file:
        file_content = csv.reader(file, delimiter="\t")
        header = next(file_content)
        header = [to_snake_case(col) for col in header]

        data = []
        for row in file_content:
            data.append(row)

    # Check if there are duplicate column names
    duplicate = [k for k, v in Counter(header).items() if v > 1]
    if len(duplicate) > 0:
        print("\nERROR: Header has duplicate columns.")
        print(f"\u2022 Duplicate columns: '{duplicate}'.")
        exit()

    outfile = path if overwrite else "_fix".join(os.path.splitext(path))
    with open(outfile, mode="w") as file:
        file_content = csv.writer(file, delimiter="\t")
        file_content.writerow(header)
        file_content.writerows(data)


def to_snake_case(string):
    """Converts a string to snake case

    Args:
        string (str): A string.
    """
    # Convert whitespace and punctuation to underscore
    regex_punc = re.compile(f"[{whitespace}{re.escape(punctuation)}]")
    string = re.sub(regex_punc, "_", string)
    
    # Add underscore before capitalized letters
    string = re.sub(r"(?<!^)(?=[A-Z])", "_", string)
    
    # Convert multiple underscores to one underscore
    string = re.sub("_+", "_", string)
    
    # Lowercase everything
    return(string.lower())

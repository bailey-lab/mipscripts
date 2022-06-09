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
    # Define whitespace and punctuation regex string
    space_punc = f"[{whitespace}{re.escape(punctuation)}]"
    
    # Remove underscore from search string
    space_punc = re.sub("_", "", space_punc)
    regex_space_punc = re.compile(space_punc)

    # Warn user that string is being converted to snake case
    if re.search(regex_space_punc, string) or re.search("[A-Z]", string):
        print(f"Converting `{string}` to snake case.")
    
    # Convert whitespace and punctuation to underscore
    string = re.sub(regex_space_punc, "_", string)

    # Add underscore before capitalized letters
    string = re.sub(r"(?<!^)(?=[A-Z])", "_", string)

    # Convert multiple underscores to one underscore
    string = re.sub("_+", "_", string)

    # Lowercase everything
    return string.lower()


def make_dir(dir, overwrite=True):
    """Checks if a directory can be made and creates it.

    Args:
        dir (str): Path to the directory.
        overwrite (bool): A boolean indicating whether to overwrite an existing
            directory.
    """
    if not isinstance(dir, str) or not dir:
        print(f"ERROR: Bad directory name of `{dir}`.")
        exit()
    elif os.path.isdir(dir) and not overwrite:
        print(f"ERROR: Directory `{dir}` already exists.")
        print("\u2022 Consider deleting the directory or renaming it.")
        exit()
    elif os.path.isdir(dir) and overwrite:
        os.system("rm -r " + dir)
        os.system("mkdir " + dir)
    else:
        os.system("mkdir " + dir)

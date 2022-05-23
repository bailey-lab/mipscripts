import csv
import os

def convert_header_snakecase(path, overwrite=False):
    """Converts header column names of a file to snakecase.

    Args:
        path (str): The path to the file.
        overwrite (bool): A boolean indicating whether to overwrite the original
            file. If `False`, a new file will be generated with the same name
            as the original file but with `_fix` appended to it.
    """
    with open(path) as file:
        file_content = csv.reader(file, delimiter="\t")
        header = next(file_content)
        header = [col.lower().replace(' ', '_') for col in header]

        data = []
        for row in file_content: data.append(row)

    outfile = path if overwrite else "_fix".join(os.path.splitext(path))
    with open(outfile, mode="w") as file:
        file_content = csv.writer(file, delimiter="\t")
        file_content.writerow(header)
        file_content.writerows(data)

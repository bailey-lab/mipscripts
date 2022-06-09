from mipscripts import utils
import csv
import os


# TODO - this should be changed to write_fastq_from_mergedsheet
def merge_sampleset_write_fastq(
    mergedsheet, newfastqdir, skipfastqwrite, skipbadfastqs
):
    """Writes FASTQ files to a directory from a merged sample sheet.

    The merged sample sheet can be created using the `merge_samplesheet` 
    function. This function allows for more complicated merges of data that may
    require intermediate manipulation of the merged sample sheet file. For
    instance, controls may be deleted or sample sets renamed. After manipulation
    of the merged sample sheet, you may run this function to combine the FASTQ
    files into one directory.

    Args:
        mergedsheet (str): A path to the merged sample sheet.
        newfastqdir (str): The path of the new FASTQ directory where files will 
            be saved.
        skipfastqwrite (bool): A boolean indicating whether to skip the writing
            of FASTQ files. This is equivalent to a dry run.
        skipbadfastqs (bool): A boolean indicating whether to skip bad FASTQs
            without throwing errors.
    """
    print(f"CREATING NEW FASTQ DIRECTORY: `{newfastqdir}`")
    utils.make_dir(newfastqdir)
    # total_bad = 0
    # total = good = 0

    with open(mergedsheet, newline="") as items_file:
        items_reader = csv.DictReader(items_file, delimiter="\t")
        # items_total = 0
        for item in items_reader:
            fastqs = item["fastq"].split(",")
            name = "{}-{}-{}".format(
                item["sample_name"], item["sample_set"], item["replicate"]
            )
            if not skipfastqwrite:
                print("...", name, "...")
            writeoperator = " > "
            if len(fastqs) % 2 == 0 and len(fastqs) > 1:
                for i in range(0, len(fastqs), 2):
                    if not skipfastqwrite:
                        os.system(
                            "cat  "
                            + fastqs[i]
                            + writeoperator
                            + newfastqdir
                            + "/"
                            + name
                            + "_R1_001.fastq.gz"
                        )
                        os.system(
                            "cat  "
                            + fastqs[i + 1]
                            + writeoperator
                            + newfastqdir
                            + "/"
                            + name
                            + "_R2_001.fastq.gz"
                        )
                    writeoperator = " >> "  # now appending after initial files
            else:
                if skipbadfastqs:
                    print("WARN bad fasta paths for ", name, "(", fastqs, ")")
                else:
                    print("ERROR: Bad FASTQ", name, fastqs)
                    print("\u2022 Skip bad FASTQs with '--skipbadfastqs'.")
                    exit(1)

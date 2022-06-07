import csv
import os


# TODO - this should be changed to write_fastq_from_mergedsheet
def merge_sampleset_write_fastq(
    mergedsheet, newfastqdir, skipfastqwrite, skipbadfastqs
):
    print("CREATING and CLEANING NEW FASTQ DIR:", newfastqdir)
    print(newfastqdir)
    if len(newfastqdir) == 1:
        os.system("mkdir " + newfastqdir)
        os.system(f"rm {newfastqdir}/*")  # too lazy to check if exist
    else:
        print(f"ERROR: bad directory name ({newfastqdir})")
        exit()
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

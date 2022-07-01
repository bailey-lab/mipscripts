from alive_progress import alive_it
from mipscripts import utils
import collections
import csv
import gzip
import os
import os.path
import pandas as pd
import re
import statistics


def seqrun_stats(samplesheet, maingrp, subgrp):
    """Computes performance statistics for a sequencing run.

    Reads the R1 FASTQ files that correspond to a sample sheet and counts the 
    number of reads for each sample. Prints summary statistics on the number
    of reads and creates a new sample sheet file that contains the original
    data and a new column `read_count` with the number of reads per sample.

    Args:
        samplesheet (list): A list of paths to sample sheets.
        maingrp (str): The primary grouping for the data.
        subgrp (str): A secondary grouping for the data.
    """

    for sheet in samplesheet:
        # Convert header to snake case
        utils.header_to_snake_case(path=sheet, overwrite=True)

        # Initialize a counter for the number of sample sets
        num_sample_sets = collections.Counter()

        # Initialize the grouping dictionary
        groupings = {}

        # Load the sample sheet header and data into memory for fast retrieval
        # Also update the number of sample per sample set
        sampledict = []
        header = []
        with open(sheet) as samplefile:
            dictreader = csv.DictReader(samplefile, delimiter="\t")
            header = dictreader.fieldnames
            for row in dictreader:
                sampledict.append(row)
                num_sample_sets.update({row["sample_set"]: 1})

        # Error if arguments are not correct
        if maingrp not in header:
            print("ERROR: '--maingrp' invalid.")
            print(f"\u2022 User entered '{maingrp}'.")
            print(f"\u2022 Available options: {header}.")
            exit()
        if subgrp and subgrp not in header:
            print("ERROR: '--subgrp' invalid.")
            print(f"\u2022 User entered '{maingrp}'.")
            print(f"\u2022 Available options: {header}.")
            exit()

        print("####################### FASTQS / READS ########################")
        fastqdir = os.path.join(os.path.dirname(sheet), "fastq")
        fastqs = os.listdir(fastqdir)
        fqss = [f for f in fastqs if "_R1_" in f]

        # Count the number of fastq reads in each fastq file
        fastqlen = {}
        bar = alive_it(fqss)
        for fq in bar:
            reads = 0
            with gzip.open(f"{fastqdir}/{fq}", mode="rt") as f:
                for line in f:
                    if line.startswith("@"):
                        reads += 1
            fastqlen[fq] = reads
            bar.title("Counting the number of FASTQ reads")

        # Save all fastq file names
        fqfiles = fastqlen.keys()

        # For each row in our data, check if there are multiple fastqs and, if
        # not, set the read count
        for row in sampledict:
            # Define identifier to compare fastqs
            fqname = "{}-{}-{}_".format(
                row["sample_name"], row["sample_set"], row["replicate"]
            )

            # Use regex to check if the fastq name appears multiple times
            regex = re.compile("^" + fqname)
            fqmatchs = list(filter(regex.match, fqfiles))

            if len(fqmatchs) > 1:
                # Throw an error
                print(
                    "ERROR: More than one FASTQ for the sampleset repetition.",
                    fqmatchs,
                    fqname,
                )
                print(
                    "\u2022 This may be because all the fastqs were not erased before"
                )
                print(
                    "  an updated sample sheet was demultiplexed (different S###)."
                )
                exit()
            elif len(fqmatchs) == 1:
                # Set read count
                row["read_count"] = fastqlen[fqmatchs[0]]

            # Create groupings
            # TODO - we don't do anything with the grouping variables
            groupings[row[maingrp]] = {}
            if subgrp:
                groupings[row[maingrp]][row[subgrp]] = 1

        # If there were undetermined reads let the user know
        if "Undetermined_S0_R1_001.fastq.gz" in fastqlen:
            cnt = fastqlen["Undetermined_S0_R1_001.fastq.gz"]
            print(f"There were {cnt} undetermined reads from demultipelxing.")

        # Create updated sample sheet
        print("Creating an updated samplesheet with read counts.")
        sheet_root = os.path.splitext(sheet)[0]
        read_cnt_path = f"{sheet_root}_readcnt.tsv"
        with open(read_cnt_path, mode="wt") as samplesheetout:
            header.append("read_count")
            dictwriter = csv.DictWriter(
                samplesheetout, delimiter="\t", fieldnames=header
            )
            dictwriter.writeheader()
            for row in sampledict:
                dictwriter.writerow(row)

        print("########################### SAMPLES ###########################")
        print("PER SAMPLE SET SUMMARY:")
        # print(f"SAMPLE GROUPINGS AND # SAMPLES: {groupings}")

        # Print some summary stats on reads for each sample set.
        # Read the read count file in
        read_cnt_data = pd.read_csv(read_cnt_path, sep="\t")
        for key in num_sample_sets:
            # Print sample set name and number of samples
            print(f"Sample set {key} ({num_sample_sets[key]:,} samples):")

            # Filter data to sample set and isolate reads
            filter_data = read_cnt_data[read_cnt_data.sample_set == key]
            num_samples = filter_data.shape[0]
            reads = filter_data.read_count

            # Print summary stats
            print(f"... TOTAL read pairs: {sum(reads):,}")
            print(f"... MEAN read pairs: {statistics.mean(reads):,}")
            print(f"... MEDIAN read pairs: {statistics.median(reads):,}")

            # Deciles
            if num_samples > 1:
                quantiles = [
                    f"{round(q):,}" for q in statistics.quantiles(reads)
                ]
                print(f"... QUANTILES for read pairs: {quantiles}")

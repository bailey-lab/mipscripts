__version__ = "v0.3.4"

import argparse
from .seqrun_stats import seqrun_stats
from .merge_sampleset_write_fastq import merge_sampleset_write_fastq
from .merge_sampleset import merge_sampleset

# Define parent argument parser and subparser
parser = argparse.ArgumentParser(
    prog="mipscripts",
    description="Utilities complementing MIPTools pipelines",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "--version",
    action="version",
    version="%(prog)s " + __version__,
)
subparsers = parser.add_subparsers(title="subcommands")

# Argument parser for seqrun_stats
parser_seqrun = subparsers.add_parser(
    "seqrun_stats",
    description="Compute performance statistics for a sequencing run",
    epilog="Note: currently this simply scans R1 FASTQs",
    help="compute performance statistics for a sequencing run",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser_seqrun.add_argument(
    "-s",
    "--samplesheet",
    action="append",
    required=True,
    help="sample sheet to run on",
)
parser_seqrun.add_argument(
    "-g",
    "--maingrp",
    default="sample_set",
    help="main grouping for stats",
)
parser_seqrun.add_argument(
    "-u",
    "--subgrp",
    help="sub grouping for stats",
)
parser_seqrun.set_defaults(func=seqrun_stats)

# Argument parser for merge_sampleset_write_fastq
parser_fastq = subparsers.add_parser(
    "merge_sampleset_write_fastq",
    description="Write FASTQ files using a merged sample sheet",
    help="write FASTQ files using a merged sample sheet",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser_fastq.add_argument(
    "-s",
    "--mergedsheet",
    required=True,
    help="merged sample sheet file paths",
)
parser_fastq.add_argument(
    "-d",
    "--newfastqdir",
    default="mergedfastq",
    help="name of the new FASTQ directory",
)
parser_fastq.add_argument(
    "-n",
    "--skipfastqwrite",
    action="store_true",
    help="skip the writing of FASTQ files. This is equivalent to a dry run.",
)
parser_fastq.add_argument(
    "-b",
    "--skipbadfastqs",
    action="store_true",
    help="skip bad FASTQs without throwing errors",
)
parser_fastq.set_defaults(func=merge_sampleset_write_fastq)

# Argument parser for merge_sampleset
parser_merge = subparsers.add_parser(
    "merge_sampleset",
    description="Flexibly merge sample sets together",
    epilog=(
        "Note: we recommend using a unique name for the outputted merged sheet"
        "based on the set/probe/and other parameters"
    ),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    help="flexibly merge sample sets together",
)
parser_merge.add_argument(
    "-S",
    "--set",
    action="append",
    required=True,
    help="sample set to aggregate",
)
parser_merge.add_argument(
    "-p",
    "--probe",
    action="append",
    required=True,
    help="probe sets to include in merge",
)
parser_merge.add_argument(
    "-s",
    "--sheet",
    action="append",
    help="sample sheet file paths (dir fastq must be in same dir as sheet)",
)
parser_merge.add_argument(
    "-o",
    "--mergeon",
    help="fields to merge on",
    default="sample_name-sample_set-replicate",
)
parser_merge.add_argument(
    "-m",
    "--newsheet",
    help="name of new merged sample sheet",
    default="mergedsheet.tsv",
)
parser_merge.add_argument(
    "-d",
    "--newfastqdir",
    help="name of new FASTQ directory",
    default="mergedfastq",
)
parser_merge.add_argument(
    "-e",
    "--exclude",
    action="append",
    help="exclude (samples) matching text pattern",
)
parser_merge.add_argument(
    "-a",
    "--addcolumn",
    action="append",
    help="add unexpected column name",
)
parser_merge.add_argument(
    "-n",
    "--skipfastqwrite",
    action="store_true",
    help="dry run skips FASTQ merge/write but makes tsv",
)
parser_merge.add_argument(
    "-r",
    "--renamereplicates",
    action="store_true",
    help="renumber the replicates based on order",
)
parser_merge.add_argument(
    "-c",
    "--collapse",
    action="store_true",
    help="collapse to unique values in columns (not needed for downstream miptools)",
)
parser_merge.add_argument(
    "-i",
    "--ignorereplicateredundancy",
    action="store_true",
    help="collapse multiple replicates within a samplesheet (error)",
)
parser_merge.set_defaults(func=merge_sampleset)

# Parse arguments and call subcommand
args = vars(parser.parse_args())
func = args.pop("func", None)
func(**args)

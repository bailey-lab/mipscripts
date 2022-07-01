from mipscripts import utils
import copy
import csv
import os

# Define global variables
global_samples_tsv_header = [
    "sample_name",
    "sample_set",
    "replicate",
    "probe_set",
    "fw",
    "rev",
    "library_prep",
]
global_samples_tsv_header_optional = [
    "capture_plate",
    "quadrant",
    "capture_plate_row",
    "capture_plate_column",
    "384_column",
    "sample_plate",
    "sample_plate_column",
    "sample_plate_row",
    "FW_plate",
    "REV_plate",
    "owner",
]
global_samples_tsv_oldnames = {
    "Library Prep": "library_prep",
    "384 Column": "384_column",
    "row": "sample_plate_row",
    "column": "sample_plate_column",
    "Sample Set": "sample_set",
}


def merge_sampleset(
    set,
    probe,
    sheet,
    mergeon,
    newsheet,
    newfastqdir,
    exclude,
    addcolumn,
    skipfastqwrite,
    renamereplicates,
    collapse,
    ignorereplicateredundancy,
):
    """Flexibly merges sample sets together into one merged sample sheet.

    Merge multiple sample sheets together, selecting the sample sets and probes
    of interest. The user may also create more advanced merges by specifying
    samples matching a pattern or selecting other fields to merge on. In some
    cases, it may be useful to edit the merged sample sheet before combining
    FASTQ files together. To do so, the user may opt to skip the creating of
    FASTQ files and later run `merged_sampleset_write_fastq`.

    Args:
        set (str): The sample sets to aggregate in the merge.
        probe (str): The probe sets to include in the merge.
        sheet (str): The path to the sample sheets that will be merged.
        mergeon (str): The fields to merge on. Fields should be separated by a
            dash.
        newsheet (str): The path to the new merged sample sheet.
        newfastqdir (str): The path to the new FASTQ directory.
        exclude (str): Exclude samples matching a text pattern.
        addcolumn (str): An unexpected column name to be added.
        skipfastqwrite (bool): Whether to skip the writing of new FASTQs. This
            is equivalent to a dry run. By setting this variable to `True`, the
            user has more flexibility in merging FASTQs. The user may use
            `merge_sampleset_write_fastq` following editing of the resulting
            sample sheet.
        renamereplicates (bool): Whether to renumber the replicates based on
            order.
        collapse (bool): Whether to collapse to unique values in columns.
        ignorereplicateredundancy (bool): Whether to collapse multiple
            replicates within a sample sheet.
    """

    # TODO: set default mergedsheet and fastq directory to set_probe
    mergeonfields = mergeon.split("-")
    if addcolumn is None:  # this should be handled by argparse but not!
        addcolumn = []

    # Convert header to snake case
    for f in sheet:
        utils.header_to_snake_case(path=f, overwrite=True)

    total_samples = []
    merged = {}
    # replicate={}
    header = global_samples_tsv_header[:]  # copy, modify
    header_optional = global_samples_tsv_header_optional[
        :
    ]  # copy, don't modify
    headers_to_rename = global_samples_tsv_oldnames  # don't modify

    print(f"CREATING NEW FASTQ DIRECTORY: `{newfastqdir}`")
    utils.make_dir(newfastqdir)
        
    print("PROCESSING SAMPLE SHEETS...")
    for thefile in sheet:
        print(f"LOADING FASTQ DIR {thefile} ...")
        print(header)

        fastqdir = os.path.join(os.path.dirname(thefile), "fastq")
        fastqs = os.listdir(fastqdir)
        fastqs = [f for f in fastqs if any(s in f for s in set)]
        # print ( "...   ",  len(fastqs)    , " fastqs in associated directory" )
        with open(thefile, newline="") as items_file:
            items_reader = csv.DictReader(items_file, delimiter="\t")
            items_total = 0
            items_kept = 0
            items_excluded = 0
            fastqs_kept = 0
            for item in items_reader:
                # rename old names to new names#
                # if items total
                for oldname in headers_to_rename:
                    if oldname in item:
                        item[headers_to_rename[oldname]] = item[oldname]
                        item.pop(oldname)
                        if items_total == 0:
                            print(
                                "... ... RENAMED old name '{}' to '{}'"
                            ).format(oldname, headers_to_rename[oldname])
                if items_total == 0:
                    # process the headerline to see what headers are present
                    # determine if we need to add another line to current header.
                    # print (item)
                    for name in item:
                        if name not in header:
                            if name in header_optional or name in addcolumn:
                                print("... ... NOTE: ", name, " added to data ")
                                header.append(name)
                            else:
                                print("ERROR: Improper column name.")
                                print(
                                    f"\u2022 '{name}' is not a proper column name."
                                )
                                print("\u2022 Use '--addcolumn' to override.")
                                exit(1)

                items_total += 1
                if not any(s == item["sample_set"] for s in set):
                    continue
                if not any(p in item["probe_set"] for p in probe):
                    continue
                # keep item with proper probeset and sampleset
                items_kept += 1

                startfastqname = "{}-{}-{}_".format(
                    item["sample_name"], item["sample_set"], item["replicate"]
                )
                item["samplelist"] = thefile
                item["replicate_old"] = item["replicate"]
                item["fastq"] = [
                    fastqdir + "/" + i
                    for i in fastqs
                    if i.startswith(startfastqname)
                ]
                item["fastq"].sort()
                fastqs_kept += len(item["fastq"])
                item["mergeon"] = []
                for f in mergeonfields:
                    item["mergeon"].append(item[f])
                item["mergeon"] = "-".join(item["mergeon"])
                if exclude and any(ele in item["mergeon"] for ele in exclude):
                    items_excluded += 1
                    continue
                # print (startfastqname)
                # find the fastqs
                total_samples.append(item)
                merged[item["mergeon"]] = 0
                # replicate[item["sample_name"]+"-"+ item["sample_set"]]=1
            print("    ", items_excluded, " excluded samples")
            print(
                "   ",
                items_kept,
                " of ",
                items_total,
                " total with ",
                fastqs_kept,
                "fastqs kept",
            )
    # print ("HEADER",header)
    # print (EXtotal_samples[0:1])
    print(
        "Total samples ",
        len(total_samples),
        " merging into ",
        len(merged),
        " samples",
    )
    total_samples.sort(key=lambda x: x["mergeon"])

    print("############################")
    print("MERGE ALL THE SAMPLE SHEET ROWS ")
    # merge the list and copy the files ####################
    replicate = {}
    merged_samples = []
    count = 0
    header = header + [
        "fastq",
        "samplelist",
        "mergeon",
        "replicate_old",
    ]  # add on additional headers to report
    for mkey in sorted(merged):
        mergeset = [i for i in total_samples if i["mergeon"] == mkey]
        count += 1
        # print (count, mkey, "records:",len(mergeset))
        mergeditem = None
        name = None

        # Create a merge record from first one and add the others
        for item in mergeset:
            if mergeditem is None:
                mergeditem = copy.deepcopy(item)
                # make sure merged item has all the potential names from header
                # set file name and replicated based on new name and sampleset#
                for f in header:
                    if f in mergeditem:
                        mergeditem[f] = [mergeditem[f]]
                    else:
                        mergeditem[f] = ["NA"]
            else:
                # print ("MERGEDITEM", mergeditem)
                for f in header:
                    if f in item:
                        mergeditem[f].append(item[f])
                    else:
                        mergeditem[f].append("NA")

        # Collapse in mergeitem sampleset, name and mergeon and error if not the same
        for f in ["sample_name", "sample_set", "mergeon"]:
            mergeditem[f] = list(set(mergeditem[f]))
            if len(mergeditem[f]) == 1:
                mergeditem[f] = mergeditem[f][0]
            else:
                print(
                    "!ERROR! ", mergeditem[f], "for", f, " is not consistent!"
                )
                exit(1)

        # Set the replicate
        name = mergeditem["sample_name"] + "-" + mergeditem["sample_set"]
        mergeditem["replicate_old"] = mergeditem["replicate"]
        mergeditem["replicate"] = mergeditem["replicate"][0]
        if renamereplicates:
            # rename replicates based on name-set
            if name in replicate:
                replicate[name] += 1
            else:
                replicate[name] = 1
            mergeditem["replicate"] = replicate[name]
        else:
            # replicate name based on name-set-replicate exists than error
            name += "-" + mergeditem["replicate"]
            if name in replicate:
                print("!ERROR!", name, ") already exists!")
                exit(1)
            else:
                replicate[name] = 1

        # Collapse probesets
        mergeditem["probe_set"] = ",".join(mergeditem["probe_set"]).split(",")
        mergeditem["probe_set"] = list(set(mergeditem["probe_set"]))
        mergeditem["probe_set"].sort()
        mergeditem["probe_set"] = ",".join(mergeditem["probe_set"])
        merged_samples.append(mergeditem)
    print("################")
    print("MERGE THE FASTQS...")
    for mergeditem in merged_samples:
        # Merge the fastqs
        # Note that this is slow
        fastqname = "{}-{}-{}".format(
            mergeditem["sample_name"],
            mergeditem["sample_set"],
            mergeditem["replicate"],
        )
        writeoperator = " > "  # initial write
        for pair in mergeditem["fastq"]:
            if len(pair) == 0:
                # print ("emptty value")
                continue
            elif (len(pair) % 2) == 0:  # 2,4,6,8
                # expectation is each replicate record has pair of fastqs (R1 & R2)
                # potential to have more than one replicate in a sequencing run but that is weird
                # print (pair)
                if len(pair) > 2:
                    if ignorereplicateredundancy:
                        print(
                            "WARN: (",
                            mergeditem["mergeon"],
                            ") just taking one pair of fastqs -- likely duplicate or bad replicate numbers!",
                        )
                    else:
                        print(mergeditem)
                        print(
                            "!ERROR! More than just a single pair for replicate!",
                            pair,
                            "\nThis could be due to not properly numbering replicates!",
                            "\nOr multiple demultiplexing runs into same directory?",
                            "\nOr because you meant to...",
                            "\nTo override use --ignorereplicateredundancy!",
                        )
                        badcount = len(
                            [m for m in total_samples if len(m["fastq"]) > 2]
                        )
                        print("NOTE: ", badcount, "have more than  2 fastqs!")
                        print(
                            "NOTE: --ignorereplicateredundancy tosses the extra fastqs"
                        )
                        exit(1)

                if not skipfastqwrite:  # this skips lengthy process of writing
                    print(
                        "...",
                        mergeditem["mergeon"],
                        "...",
                        pair[0][-65:],
                        "...",
                    )
                    os.system(
                        "cat  "
                        + pair[0]
                        + writeoperator
                        + newfastqdir
                        + "/"
                        + fastqname
                        + "_R1_001.fastq.gz"
                    )
                    os.system(
                        "cat  "
                        + pair[1]
                        + writeoperator
                        + newfastqdir
                        + "/"
                        + fastqname
                        + "_R2_001.fastq.gz"
                    )
                    writeoperator = " >> "  # now appending after initial files
                # else:
                #  for i in range (0, len(pair),2):
                #     print (pair[i],"\n   ", pair[i+1])
            else:
                print("ODD PAIR ERROR ", pair)
                print(" could a R1 or R2 fastq file been deleted???")
                exit(1)
            writeoperator = " >> "  # now appending
        # Make fields non-redundant (always make probeset)
        mergeditem["fastq"] = [
            fq for sublist in mergeditem["fastq"] for fq in sublist
        ]
        # print (mergeditem)
    # if collapse :
    #   uniquefields=fieldstocollapse
    #  for f in uniquefields :
    #  if f in mergeditem:
    #      #print (f, mergeditem[f])
    #     mergeditem[f]= list(set(mergeditem[f]))
    # for f in fieldstocollapse :
    #  if f in mergeditem:
    #    mergeditem[f]= ",".join(mergeditem[f])

    # print ("NEW COPY ", mergeditem)
    # write to a file ################################
    # write header and additional data from merge
    print("WRITING TSV: ", newsheet)
    with open(newsheet, "wt") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(header)
        for item in merged_samples:
            row = []
            for f in header:
                if type(item[f]) is list:
                    if collapse:
                        item[f] = list(set(item[f]))
                    item[f] = ",".join(item[f])
                row.append(item[f])
            tsv_writer.writerow(row)

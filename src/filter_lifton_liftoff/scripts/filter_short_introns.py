#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to filter GFF3 on intron size

"""

# import libraries
import argparse
import os
import re
import logging
import sys
from collections import defaultdict
import statistics

# get script name
script = os.path.basename(sys.argv[0])

# get the GFF3 attributes
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTE = range(9)

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)


def get_id(file, line, attribute, field):
    """
    Extract query from GFF3 attribute column
    """
    pattern = field + "=([^;]+)"
    try:
        # Check for GFF3 file
        id_search = re.search(pattern, attribute)
        id_field = None
        if id_search:
            id_field = id_search.group(1)
        else:
            logging.error(
                f"Error: Cannot extract attribute type '{field}=' from the file '{file}'. Please check if the attribute is present in line:'\n{line}\n"
            )
    except AttributeError as err:
        logging.error(
            f"Error: {err}. '{field}' cannot be extracted from the file '{file}' line below\n{line}"
        )

    try:
        id_field
    except NameError as err:
        logging.error(f"Error: Cannot find {field} from the file '{file}', exiting..")
        sys.exit(1)

    return id_field


def get_gff3_info(filename):
    logging.info("Processing file {filename} ...".format(filename=filename))
    gene_info = defaultdict(dict)
    mrna_info = defaultdict(dict)
    exon_info = defaultdict(list)
    with open(filename, "r") as filehandle:
        for line in filehandle:
            line = line.rstrip("\n")
            if re.match(r"^\s*$", line) or line.startswith("#"):
                # skip directives and comments
                continue
            x = line.split("\t")
            if len(x) != 9:
                logging.error(
                    "Error: Not a standard GFF3 9 column line, see below:\n{line}. Exiting...".format(
                        line=line
                    )
                )
                sys.exit(1)

            start, end = (
                (int(x[START]), int(x[END]))
                if int(x[START]) < int(x[END])
                else (int(x[END]), int(x[START]))
            )

            if x[TYPE] in ["gene"]:
                gene_id = get_id(filename, line, x[ATTRIBUTE], "ID")
                if gene_id in gene_info:
                    logging.error(
                        "Error: Duplicate gene encountered '{gene_id}'. Please make sure that the GFF3 file is sorted using genometools. Exiting...".format(
                            gene_id=gene_id
                        )
                    )
                    sys.exit(1)
                else:
                    gene_info[gene_id]["start"] = start
                    gene_info[gene_id]["end"] = end
            elif x[TYPE] in ["mRNA"]:
                mrna_id = get_id(filename, line, x[ATTRIBUTE], "ID")
                parent_id = get_id(filename, line, x[ATTRIBUTE], "Parent")
                if mrna_id in mrna_info:
                    logging.error(
                        "Error: Duplicate mRNA encountered '{mrna_id}'. Please make sure that the GFF3 file is sorted using genometools. Exiting...".format(
                            mrna_id=mrna_id
                        )
                    )
                    sys.exit(1)
                else:
                    mrna_info[mrna_id]["parent_id"] = parent_id
                    mrna_info[mrna_id]["start"] = start
                    mrna_info[mrna_id]["end"] = end

            elif x[TYPE] in ["exon"]:
                parent_id = get_id(filename, line, x[ATTRIBUTE], "Parent")
                exon_info[parent_id].append([start, end])
    return gene_info, mrna_info, exon_info


def remove_duplicate_and_merge(exon_info):
    logging.info("Remove duplicate and merge overlapping exons ...")
    # remove duplicate exon coordinates
    for transcript in exon_info:
        exon_info[transcript] = [
            list(x) for x in set(tuple(x) for x in exon_info[transcript])
        ]
    # merge overlapping exons
    for transcript in exon_info:
        exon_info[transcript].sort(key=lambda interval: interval[0])
        merged = []
        for current in exon_info[transcript]:
            if not merged:
                merged.append(current)
            else:
                previous = merged[-1]
                if current[0] <= previous[1]:  # Overlap
                    previous[1] = max(previous[1], current[1])  # Merge
                else:
                    merged.append(current)
        exon_info[transcript] = merged
    return exon_info


def compute_intron_stats(exon_info):
    logging.info("Compute intron sizes ...")
    intron_info = defaultdict(list)
    for transcript, coords in exon_info.items():
        coords.sort(key=lambda interval: interval[0])
        # print(transcript, coords)
        if len(coords) == 1:  # For mono exonic models
            intron_info[transcript] = []
        else:
            for i, _ in enumerate(coords[:-1]):
                intron_info[transcript].append([coords[i][1] + 1, coords[i + 1][0] - 1])
                # print(f"i:{i} {coords[i][1] + 1 } {coords[i+1][0] - 1}")

    # logging.info(
    #     "# INFO: Compute max intron size and intron sizes at ends ...")
    intron_stats_info = defaultdict(dict)
    for transcript, coords in intron_info.items():
        coords.sort(key=lambda interval: interval[0])
        if len(coords) == 0:  # for mono exonic models
            intron_stats_info[transcript]["start"] = 0
            intron_stats_info[transcript]["end"] = 0
            intron_stats_info[transcript]["min"] = 0
            intron_stats_info[transcript]["max"] = 0
        else:
            intron_stats_info[transcript]["start"] = coords[0][1] - coords[0][0] + 1
            intron_stats_info[transcript]["end"] = coords[-1][1] - coords[-1][0] + 1
            intron_stats_info[transcript]["min"] = min(
                [coord[1] - coord[0] + 1 for coord in coords]
            )
            intron_stats_info[transcript]["max"] = max(
                [coord[1] - coord[0] + 1 for coord in coords]
            )
    # Get intron stats
    min_intron = min(
        value
        for k, v in intron_stats_info.items()
        for x, value in v.items()
        if x in "min"
    )
    max_intron = max(
        value
        for k, v in intron_stats_info.items()
        for x, value in v.items()
        if x in "max"
    )
    mean_intron = statistics.mean(
        (
            value
            for k, v in intron_stats_info.items()
            for x, value in v.items()
            if x in "max"
        )
    )
    median_intron = statistics.median(
        (
            value
            for k, v in intron_stats_info.items()
            for x, value in v.items()
            if x in "max"
        )
    )
    min_ends_intron = min(
        value
        for k, v in intron_stats_info.items()
        for x, value in v.items()
        if x in ["start", "end"]
    )
    max_ends_intron = max(
        value
        for k, v in intron_stats_info.items()
        for x, value in v.items()
        if x in ["start", "end"]
    )
    logging.info(f"# - Minimum intron size:\t{min_intron}")
    logging.info(f"# - Maximum intron size:\t{max_intron}")
    logging.info(f"# - Mean intron size:\t{mean_intron}")
    logging.info(f"# - Median intron size:\t{median_intron}")
    logging.info(f"# - Minimum terminal introns size:\t{min_ends_intron}")
    logging.info(f"# - Maximum terminal introns size:\t{max_ends_intron}")
    return intron_info, intron_stats_info


def filter_intron(
    gene_info,
    mrna_info,
    exon_info,
    intron_stats_info,
    intron_size,
    intron_size_ends,
    detailed_discard_list_file,
):
    # logging.info("Filtering intron ...")
    discard_list = {}
    for transcript in intron_stats_info:
        # only remove multi-exonic models
        if len(exon_info[transcript]) < 2:
            continue
        if intron_size and intron_stats_info[transcript]["min"] < intron_size:
            if transcript not in discard_list:
                discard_list[transcript] = mrna_info[transcript]["parent_id"]
        if intron_size_ends and (
            intron_stats_info[transcript]["start"] < intron_size_ends
            or intron_stats_info[transcript]["end"] < intron_size_ends
        ):
            if transcript not in discard_list:
                discard_list[transcript] = mrna_info[transcript]["parent_id"]
    # print(discard_list)

    logging.info(
        "Identified {count} mRNAs (of total {mrna_count} mRNAs / {gene_count} genes) to be removed".format(
            count=len(discard_list.keys()),
            mrna_count=len(mrna_info.keys()),
            gene_count=len(gene_info.keys()),
        )
    )
    if intron_size:
        logging.info(
            "# - based on intron size <= {intron_size} bp".format(
                intron_size=intron_size
            )
        )
    if intron_size_ends:
        logging.info(
            "# - based on intron size ends <= {intron_size_ends} bp".format(
                intron_size_ends=intron_size_ends
            )
        )

    # logging.info("INFO: Removing models that failed intron size, if any ...")
    # add detailed discard list
    # mrna id , gene id, exon coords, min intron size, max intron size, exon lines
    with open(detailed_discard_list_file, "w") as f:
        f.write(
            "\t".join(
                [
                    "mrna_id",
                    "gene_id",
                    "exon_coords",
                    "min_intron_size",
                    "max_intron_size",
                ]
            )
            + "\n"
        )
        for discard_models in discard_list.keys():
            f.write(
                "\t".join(
                    [
                        discard_models,
                        discard_list[discard_models],
                        str(exon_info[discard_models]),
                        str(intron_stats_info[discard_models]["min"]),
                        str(intron_stats_info[discard_models]["max"]),
                    ]
                )
                + "\n"
            )

    filtered_mrna_info = mrna_info.copy()
    for mrna in discard_list:
        if mrna in filtered_mrna_info:
            del filtered_mrna_info[mrna]

    logging.info("Filter failed intron models and readjust gene coordinates ...")
    updated_gene_info = defaultdict(dict)
    for mrna in filtered_mrna_info:
        start = filtered_mrna_info[mrna]["start"]
        end = filtered_mrna_info[mrna]["end"]
        parent_id = filtered_mrna_info[mrna]["parent_id"]
        if parent_id not in updated_gene_info:
            updated_gene_info[parent_id]["start"] = start
            updated_gene_info[parent_id]["end"] = end
        else:
            prev_start = updated_gene_info[parent_id]["start"]
            prev_end = updated_gene_info[parent_id]["end"]
            if int(start) < int(prev_start):
                updated_gene_info[parent_id]["start"] = start
            if int(end) > int(prev_end):
                updated_gene_info[parent_id]["end"] = end
    return discard_list, filtered_mrna_info, updated_gene_info


def process_output(
    filename,
    source,
    intron_info,
    add_intron,
    gene_info,
    mrna_info,
    discard_list,
    discard_list_mapping_file,
    detailed_discard_list_file,
    filtered_mrna_info,
    updated_gene_info,
    output_gff,
):
    # write output GFF3
    with open(output_gff, "w") as f:
        f.write("##gff-version 3\n")
        skip_first = False
        with open(filename, "r") as filehandle:
            for line in filehandle:
                line = line.rstrip("\n")
                if re.match(r"^\s*$", line) or line.startswith("#"):
                    # skip directives and comments
                    pass
                else:
                    x = line.split("\t")
                    # change source if requested
                    if source:
                        x[SOURCE] = source
                    if x[TYPE] in ["gene"]:
                        gene_id = get_id(filename, line, x[ATTRIBUTE], "ID")
                        if gene_id in updated_gene_info:
                            if skip_first:
                                f.write("###\n")
                            skip_first = True
                            new_start = updated_gene_info[gene_id]["start"]
                            new_end = updated_gene_info[gene_id]["end"]
                            fline = [
                                x[SEQID],
                                x[SOURCE],
                                x[TYPE],
                                str(new_start),
                                str(new_end),
                                x[SCORE],
                                x[STRAND],
                                x[PHASE],
                                x[ATTRIBUTE],
                            ]
                            f.write("\t".join(fline) + "\n")

                    elif x[TYPE] in ["mRNA"]:
                        mrna_id = get_id(filename, line, x[ATTRIBUTE], "ID")
                        if mrna_id in filtered_mrna_info:
                            f.write("\t".join(x) + "\n")
                            # add intron GFF3 lines, if true
                            if add_intron:
                                if mrna_id in intron_info:
                                    if not len(intron_info[mrna_id]) == 0:
                                        for i, coord in enumerate(intron_info[mrna_id]):
                                            intronLine = [
                                                x[SEQID],
                                                x[SOURCE],
                                                "intron",
                                                str(coord[0]),
                                                str(coord[1]),
                                                ".",
                                                x[STRAND],
                                                ".",
                                                f"ID={mrna_id}.intron{i+1};Parent={mrna_id}",
                                            ]
                                            f.write("\t".join(intronLine) + "\n")
                    elif x[TYPE] in [
                        "exon",
                        "CDS",
                        "five_prime_UTR",
                        "three_prime_UTR",
                        "five_prime_utr",
                        "three_prime_utr",
                    ]:
                        parent_id = get_id(filename, line, x[ATTRIBUTE], "Parent")
                        if parent_id in filtered_mrna_info:
                            f.write("\t".join(x) + "\n")
        f.write("###\n")

    # write discard list
    with open(discard_list_mapping_file, "w") as f:
        f.write("mrna_id\tgene_id\n")
        for discard_models in discard_list.keys():
            # f.write(discard_models + "\n")
            f.write(f"{discard_models}\t{discard_list[discard_models]}\n")

    logging.info("SUMMARY")
    logging.info(f"# Input file:\t{filename}")
    logging.info(f"# Input genes:\t{len(gene_info.keys())}")
    logging.info(f"# Input mRNAs:\t{len(mrna_info.keys())}")
    logging.info(f"# mRNAs removed based on intron filter:\t{len(discard_list.keys())}")
    logging.info(
        f"# genes removed based on intron filter:\t{len(set(discard_list.values()))}"
    )
    logging.info(
        f"# Output genes:\t{len(updated_gene_info.keys())} [{(100 * len(updated_gene_info.keys()) / len(gene_info.keys())):.2f}%] "
    )
    logging.info(
        f"# Output mRNAs:\t{len(filtered_mrna_info.keys())} [{(100 * len(filtered_mrna_info.keys()) / len(mrna_info.keys())):.2f}%]"
    )

    with open(detailed_discard_list_file, "a") as f:
        f.write(f"\n\n")
        f.write(f"Category,SubCategory,Count,Unique TID Count\n")
        f.write(
            f"discarding,Short intron,{len(discard_list.keys())},{len(discard_list.keys())}\n"
        )
        f.write(f"Checks complete.\n")

    logging.info(
        f"# Discarded models mapping file written to file: {discard_list_mapping_file}"
    )
    logging.info(
        f"# Discarded models detailed information written to file: {detailed_discard_list_file}"
    )
    logging.info(f"# Filtered output GFF3 written to file: {output_gff}")


def main():

    parser = argparse.ArgumentParser(
        description="Script to filter GFF3 on intron size for GFF3 structure: [gene]->[mRNA]->[exon|CDS|five_prime_UTR|three_prime_UTR|five_prime_utr|three_prime_utr]"
    )
    parser.add_argument("file", help="Provide GFF3 file [required]")
    parser.add_argument(
        "--source",
        default=None,
        help="Provide new source for GFF3 (default: %(default)s)",
    )
    parser.add_argument(
        "--intron_size",
        default=0,
        type=int,
        help="Provide intron size. Any models with intron size under this value will be removed (default: < %(default)s)",
    )
    parser.add_argument(
        "--intron_size_ends",
        default=0,
        type=int,
        help="Provide intron size for ends. Any models with terminal intron size under this value will be removed (default: < %(default)s)",
    )
    parser.add_argument(
        "--add_intron",
        action="store_true",
        help="Add intron type to the output GFF3 (default: %(default)s)",
    )
    # discard list file
    parser.add_argument(
        "--discard_list_mapping_file",
        required=True,
        help="Output TSV mapping mrna_id to gene_id for qualifying introns [required] (default: %(default)s)",
    )
    # detailed dicard list file
    parser.add_argument(
        "--detailed_discard_list_file",
        required=True,
        help="Output TSV listing introns meeting the filter criteria [required] (default: %(default)s)",
    )
    # output file
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output GFF3 file [required] (default: %(default)s)",
    )

    args = parser.parse_args()
    filename = args.file
    source = args.source
    intron_size = args.intron_size
    intron_size_ends = args.intron_size_ends
    add_intron = args.add_intron
    discard_list_mapping_file = os.path.abspath(args.discard_list_mapping_file)
    detailed_discard_list_file = os.path.abspath(args.detailed_discard_list_file)
    output_gff = os.path.abspath(args.output)

    # create parent dir if not exists
    os.makedirs(os.path.dirname(discard_list_mapping_file), exist_ok=True)
    os.makedirs(os.path.dirname(detailed_discard_list_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_gff), exist_ok=True)

    gene_info, mrna_info, exon_info = get_gff3_info(filename)
    exon_info = remove_duplicate_and_merge(exon_info)
    intron_info, intron_stats_info = compute_intron_stats(exon_info)
    discard_list, filtered_mrna_info, updated_gene_info = filter_intron(
        gene_info,
        mrna_info,
        exon_info,
        intron_stats_info,
        intron_size,
        intron_size_ends,
        detailed_discard_list_file,
    )
    process_output(
        filename,
        source,
        intron_info,
        add_intron,
        gene_info,
        mrna_info,
        discard_list,
        discard_list_mapping_file,
        detailed_discard_list_file,
        filtered_mrna_info,
        updated_gene_info,
        output_gff,
    )


if __name__ == "__main__":
    main()

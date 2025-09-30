#!/usr/bin/env python3

import argparse
import re
import csv


def parse_list_file(list_file):
    """
    Parse the list file to get a set of mRNA IDs.
    """
    with open(list_file, "r") as file:
        return set(line.strip() for line in file if line.strip())


def extract_parent_ids(gff_file, mrna_ids, output_file):
    """
    Extract Parent IDs for mRNA features matching the IDs in the list file.
    """
    output_data = []

    with open(gff_file, "r") as file:
        for line in file:
            # Skip header or comment lines
            if line.startswith("#"):
                continue

            # Split the GFF line into columns
            columns = line.strip().split("\t")
            if len(columns) < 9:
                continue

            feature_type = columns[2]
            attributes = columns[8]

            # Process only mRNA features
            if feature_type == "mRNA":
                # Extract the ID and Parent attributes
                id_match = re.search(r"ID=([^;]+)", attributes)
                parent_match = re.search(r"Parent=([^;]+)", attributes)

                if id_match and parent_match:
                    mrna_id = id_match.group(1)
                    parent_id = parent_match.group(1)

                    # Check if the mRNA ID matches any in the list
                    if mrna_id in mrna_ids:
                        output_data.append((mrna_id, parent_id))

    # Write the output to a TSV file without a header
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerows(output_data)

    print(f"Extraction complete. Output written to {output_file}")


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Extract Parent IDs for matching mRNA features from a GFF file."
    )
    parser.add_argument(
        "list_file", type=str, help="Path to the list file containing mRNA IDs."
    )
    parser.add_argument("gff_file", type=str, help="Path to the GFF file.")
    parser.add_argument("output_file", type=str, help="Path to the output TSV file.")
    args = parser.parse_args()

    # Parse the list file
    mrna_ids = parse_list_file(args.list_file)

    # Extract Parent IDs and write to output
    extract_parent_ids(args.gff_file, mrna_ids, args.output_file)


if __name__ == "__main__":
    main()

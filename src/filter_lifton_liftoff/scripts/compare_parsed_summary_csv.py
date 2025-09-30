#!/usr/bin/env python3

import argparse
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet
import warnings

# Suppress FutureWarning
warnings.simplefilter(action="ignore", category=FutureWarning)


def parse_gff(input_files, output_dir, features):
    for gff_file in input_files:
        output_file = os.path.join(
            output_dir, f"{os.path.basename(os.path.splitext(gff_file)[0])}.ids.tsv"
        )
        with open(gff_file, "r") as infile, open(output_file, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")
                if len(columns) < 9:
                    continue

                feature_type = columns[2]
                if feature_type not in features:
                    continue

                attributes = columns[8]
                attrs = {}
                for raw in attributes.split(";"):
                    if "=" not in raw:
                        continue
                    key, val = raw.strip().split("=", 1)
                    attrs[key] = val
                feature_id = attrs.get("ID", "")
                parent_id = attrs.get("Parent", "")

                if feature_id:
                    outfile.write(f"{feature_id}\t{parent_id}\n")


def parse_summary_csv(csv_files, output_dir, exclude_types, rm_prefix):
    for csv_file in csv_files:
        output_file = os.path.join(
            output_dir,
            f"{os.path.basename(os.path.splitext(csv_file)[0])}.ids.rejected.txt",
        )
        unique_ids = set()

        with open(csv_file, "r") as infile:
            reader = csv.reader(infile)
            next(reader)  # Skip the header
            for row in reader:
                if len(row) < 3 or not row[2].strip():
                    continue

                exclude_type = row[1].strip()
                if exclude_type in exclude_types:
                    continue

                transcript_id = row[2].strip()
                if rm_prefix and transcript_id.startswith(rm_prefix):
                    transcript_id = transcript_id[len(rm_prefix) :]
                unique_ids.add(transcript_id)

        with open(output_file, "w") as outfile:
            for transcript_id in sorted(unique_ids):
                outfile.write(f"{transcript_id}\n")


def generate_retained_lists(gff_ids_files, rejected_ids_files, output_dir):
    for gff_ids_file, rejected_ids_file in zip(gff_ids_files, rejected_ids_files):
        retained_output_file = os.path.join(
            output_dir,
            f"{os.path.basename(os.path.splitext(os.path.splitext(gff_ids_file)[0])[0])}.ids.retained.txt",
        )
        gff_ids = set()
        rejected_ids = set()

        with open(gff_ids_file, "r") as infile:
            for line in infile:
                gff_id = line.strip().split("\t")[0]
                gff_ids.add(gff_id)

        with open(rejected_ids_file, "r") as infile:
            for line in infile:
                rejected_ids.add(line.strip())

        retained_ids = gff_ids.difference(rejected_ids)

        with open(retained_output_file, "w") as outfile:
            for retained_id in sorted(retained_ids):
                outfile.write(f"{retained_id}\n")


def generate_upset_plot(
    rejected_ids_files,
    retained_ids_files,
    labels,
    base_output_name,
    output_dir,
    plot_title,
):
    ids_1 = set()
    ids_2 = set()
    retained_1 = set()
    retained_2 = set()

    with open(rejected_ids_files[0], "r") as file:
        for line in file:
            ids_1.add(line.strip())

    with open(rejected_ids_files[1], "r") as file:
        for line in file:
            ids_2.add(line.strip())

    with open(retained_ids_files[0], "r") as file:
        for line in file:
            retained_1.add(line.strip())

    with open(retained_ids_files[1], "r") as file:
        for line in file:
            retained_2.add(line.strip())

    all_ids = sorted(ids_1.union(ids_2).union(retained_1).union(retained_2))

    label_1, label_2 = labels
    data_frame = pd.DataFrame(
        {
            f"Rejected_{label_1}": [id_ in ids_1 for id_ in all_ids],
            f"Rejected_{label_2}": [id_ in ids_2 for id_ in all_ids],
            f"Retained_{label_1}": [id_ in retained_1 for id_ in all_ids],
            f"Retained_{label_2}": [id_ in retained_2 for id_ in all_ids],
        },
        index=all_ids,
    )

    presence_absence_output = os.path.join(
        output_dir, f"{base_output_name}_presence_absence_data.tsv"
    )
    data_frame.to_csv(presence_absence_output, sep="\t")
    print(f"\nPresence/Absence DataFrame saved to '{presence_absence_output}'")

    multiindex_frame = data_frame.reset_index()
    multiindex_frame = multiindex_frame.melt(
        id_vars="index", var_name="Set", value_name="Presence"
    )
    multiindex_frame = multiindex_frame[multiindex_frame["Presence"] == True]
    multiindex_frame = multiindex_frame.pivot_table(
        index="index", columns="Set", aggfunc=lambda x: 1, fill_value=0
    )
    multiindex_frame.columns = [
        col[1] if isinstance(col, tuple) else col for col in multiindex_frame.columns
    ]
    multiindex_frame.index = pd.MultiIndex.from_frame(
        multiindex_frame.reset_index(drop=True)
    )
    multiindex_frame = multiindex_frame.astype(bool)

    try:
        upset = UpSet(
            multiindex_frame, subset_size="count", intersection_plot_elements=5
        )
        fig = upset.plot()

        ax = plt.gca()
        for patch in ax.patches:
            if patch.get_height() > 0:
                ax.text(
                    patch.get_x() + patch.get_width() / 2,
                    patch.get_height() + 0.5,
                    f"{int(patch.get_height())}",
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    color="black",
                )

        if plot_title:
            plt.suptitle(plot_title, fontsize=8)

        plot_output = os.path.join(output_dir, f"{base_output_name}_upset_plot.png")
        plt.savefig(plot_output, dpi=300)
        print(f"\nUpSet plot saved to '{plot_output}'\n")
    except Exception as e:
        print("\nError during UpSet plot generation:")
        print(e)


def main():
    parser = argparse.ArgumentParser(
        description="Parse GFF and CSV files to extract and compare transcript IDs."
    )
    parser.add_argument(
        "--gff_files",
        metavar="GFF",
        type=str,
        nargs=2,
        help="Paths to two GFF files to be processed.",
    )
    parser.add_argument(
        "--csv_files",
        metavar="CSV",
        type=str,
        nargs=2,
        help="Paths to two parsed summary CSV files to be processed.",
    )
    parser.add_argument(
        "--gff_file", type=str, help="Single GFF file to process in single mode"
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        help="Single parsed summary CSV file to process in single mode",
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Enable single mode: process one GFF and one CSV file without comparisons or plots",
    )
    parser.add_argument(
        "--features",
        type=str,
        default="mRNA,transcript",
        help="Comma-separated list of features to extract (default: mRNA,transcript)",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        default="",
        help="Comma-separated list of types to exclude",
    )
    parser.add_argument(
        "--rm_prefix", type=str, default="", help="Prefix to remove from transcript IDs"
    )
    parser.add_argument(
        "--labels",
        type=str,
        default="1,2",
        help='Comma-separated list of labels for the _1 and _2 sets (default: "1,2")',
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="output",
        help='Base output name for plot and TSV file (default: "output")',
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Provide output directory (default: current directory)",
    )
    parser.add_argument(
        "--plot_title",
        type=str,
        default="",
        help="Title for the UpSet plot (default: empty)",
    )
    args = parser.parse_args()

    features = args.features.split(",")
    exclude_types = (
        [etype.strip() for etype in args.exclude.split(",")] if args.exclude else []
    )
    rm_prefix = args.rm_prefix if args.rm_prefix else None
    labels = args.labels.split(",")

    if args.single:
        if not args.gff_file or not args.csv_file:
            raise ValueError(
                "Single mode requires --gff_file and --csv_file arguments."
            )

        parse_gff([args.gff_file], args.output, features)
        parse_summary_csv([args.csv_file], args.output, exclude_types, rm_prefix)

        gff_ids_file = os.path.join(
            args.output,
            f"{os.path.basename(os.path.splitext(args.gff_file)[0])}.ids.tsv",
        )
        rejected_ids_file = os.path.join(
            args.output,
            f"{os.path.basename(os.path.splitext(args.csv_file)[0])}.ids.rejected.txt",
        )
        generate_retained_lists([gff_ids_file], [rejected_ids_file], args.output)
    else:
        if not args.gff_files or not args.csv_files:
            raise ValueError(
                "Dual mode requires --gff_files and --csv_files (two of each)."
            )

        parse_gff(args.gff_files, args.output, features)
        parse_summary_csv(args.csv_files, args.output, exclude_types, rm_prefix)

        gff_ids_files = [
            os.path.join(
                args.output,
                f"{os.path.basename(os.path.splitext(gff_file)[0])}.ids.tsv",
            )
            for gff_file in args.gff_files
        ]
        rejected_ids_files = [
            os.path.join(
                args.output,
                f"{os.path.basename(os.path.splitext(csv_file)[0])}.ids.rejected.txt",
            )
            for csv_file in args.csv_files
        ]
        generate_retained_lists(gff_ids_files, rejected_ids_files, args.output)

        retained_ids_files = [
            os.path.join(
                args.output,
                f"{os.path.basename(os.path.splitext(gff_file)[0])}.ids.retained.txt",
            )
            for gff_file in args.gff_files
        ]
        generate_upset_plot(
            rejected_ids_files,
            retained_ids_files,
            labels,
            args.output_prefix,
            args.output,
            args.plot_title,
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Generate stacked bar plots from summary CSV files."
    )
    parser.add_argument(
        "summary_csvs",
        metavar="summary_csvs",
        type=str,
        nargs="+",
        help="One or more *_summary_stats.csv files to include in the plot.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Transcript Filtering Summary",
        help="Title for the plot.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="summary_plot",
        help="Base name for output PNG plot (no extension).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output directory.",
    )
    args = parser.parse_args()

    # Color palette (sorted alphabetically by key)
    color_palette = {
        "5'UTR present with a truncated ORF": "#cc6666",
        "Assertion failure": "#c5b0d5",
        "Assertion failure start must be less than end": "#ff9896",
        "Both UTR present with truncated ORF": "#d62728",
        "CDS which straddles 2 different exons": "#8c564b",
        "Cannot reverse strand of coding transcript": "#e377c2",
        "Debords its exon": "#bcbddc",
        "Defined UTRs but no CDS feature": "#9467bd",
        "Duplicate parent feature ID": "#a65628",
        "General": "#5A5A5A",
        "Incorrect fusions of splice junctions": "#2ca02c",
        "Internal stop codons found": "#ff7f0e",
        "Invalid CDS length": "#98df8a",
        "Invalid number of coding exons": "#bcbd22",
        "Invalid start and stop of the ORF": "#17becf",
        "Overlapping CDS": "#9edae5",
        "Overlapping exons found": "#7f7f7f",
        "Redundant": "#c49c94",
        "Seqid mismatch": "#1b9e77",
        "Short intron": "#e6ab02",
        "Size under minimum": "#1f77b4",
        "Strand conflict child": "#d95f02",
        "Strand conflict gene": "#7570b3",
        "Strand conflict gene-child": "#66c2a5",
    }

    # Prepare the summary DataFrame
    summary_data = defaultdict(
        lambda: defaultdict(int)
    )  # sample -> subcategory -> count

    for csv_path in args.summary_csvs:
        if not os.path.isfile(csv_path):
            print(f"Warning: File '{csv_path}' not found, skipping.")
            continue

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"Failed to read '{csv_path}': {e}")
            continue

        # Derive sample name from filename
        # Remove "_summary_stats.csv" or "mikado_prepare_" if present from the file name for labeling
        sample = (
            os.path.basename(csv_path)
            .replace("_summary_stats.csv", "")
            .replace("mikado_prepare_", "")
        )

        for _, row in df.iterrows():
            subcat = row.get("SubCategory", "Unknown")
            count = row.get("Count", 0)
            cat = row.get("Category", "")

            if cat == "Total":
                continue  # skip the total line

            try:
                summary_data[sample][subcat] += int(count)
            except Exception:
                # If count isn't cleanly castable, skip that row with a warning.
                print(
                    f"Warning: Non-integer count '{count}' in '{csv_path}' for '{subcat}', skipping."
                )

    # Convert to DataFrame
    plot_df = pd.DataFrame(summary_data).fillna(0)

    if plot_df.empty:
        print("No valid data to plot. Exiting.")
        sys.exit(0)

    # First: alphabetical order of palette-defined subcategories that actually appear
    ordered_subcats = [k for k in sorted(color_palette.keys()) if k in plot_df.index]
    # Second: any subcategories not in the palette (unexpected/new), append alphabetically at the end
    other_subcats = sorted([idx for idx in plot_df.index if idx not in color_palette])
    ordered_index = ordered_subcats + other_subcats

    plot_df = plot_df.reindex(ordered_index)

    # Plot setup
    labels = plot_df.columns.tolist()
    fig, ax = plt.subplots(figsize=(12, 8))
    bottoms = [0] * len(labels)

    # Stack the bars in the enforced order
    for subcategory in plot_df.index:
        values = plot_df.loc[subcategory].tolist()
        ax.bar(
            labels,
            values,
            bottom=bottoms,
            label=subcategory,
            color=color_palette.get(subcategory, "#333333"),
        )
        bottoms = [i + j for i, j in zip(bottoms, values)]

    # Final plot adjustments
    plt.ylabel("Transcript Count")
    plt.xlabel("Sample")
    plt.title(args.title, loc="center", pad=20)
    plt.xticks(rotation=45, ha="right", fontsize=8)

    # Place legend below plot, ensuring it doesn't overlap x-axis
    legend = plt.legend(
        title="SubCategory",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.65),  # push it down further
        ncol=4,
        fontsize=8,
        title_fontsize=10,
        frameon=False,
    )

    # Increase space at the bottom to fit legend and labels
    plt.subplots_adjust(bottom=0.35)
    plt.tight_layout()

    # Save plot
    out_dir = os.path.abspath(args.output)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{args.output_prefix}.png")
    plt.savefig(out_path, bbox_inches="tight")  # ensure legend is included
    print(f"Plot saved to '{out_path}'")


if __name__ == "__main__":
    main()

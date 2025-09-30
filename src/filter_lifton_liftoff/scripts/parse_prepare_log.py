#!/usr/bin/env python3

import re
import csv
import sys
import os
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="Parse log files and generate summary statistics and plots."
    )
    parser.add_argument(
        "log_files",
        metavar="log_files",
        type=str,
        nargs="+",
        help="Paths to log files to be parsed.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Count Of Rejected Transcripts",
        help="Title for the generated plot.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="summary_stacked_bar_plot",
        help="Base name for the output plot file (without extension).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output directory.",
    )
    args = parser.parse_args()

    log_file_paths = args.log_files
    plot_title = args.title
    plot_output_base = os.path.join(os.path.abspath(args.output), args.output_prefix)

    # Define possible categories and subcategories
    all_combinations = [
        ("discarding", "Size under minimum"),
        ("discarding", "Internal stop codons found"),
        ("discarding", "Incorrect fusions of splice junctions"),
        ("discarding", "Both UTR present with truncated ORF"),
        ("discarding", "Defined UTRs but no CDS feature"),
        ("discarding", "CDS which straddles 2 different exons"),
        ("discarding", "Cannot reverse strand of coding transcript"),
        ("discarding", "Overlapping exons found"),
        ("discarding", "Invalid number of coding exons"),
        ("discarding", "Invalid start and stop of the ORF"),
        ("discarding", "Debords its exon"),
        ("discarding", "5'UTR present with a truncated ORF"),
        ("discarding", "General"),
        ("validation", "Assertion failure start must be less than end"),
        ("validation", "Assertion failure"),
        ("validation", "General"),
        ("excluding", "Redundant"),
        ("excluding", "General"),
    ]

    # Main patterns for the primary categories
    main_patterns = {
        "discarding": re.compile(r"\b[Dd]iscard(?:ing|ed|s)\b"),
        "excluding": re.compile(r"\b[Ee]xcluding\b"),
        "validation": re.compile(r"\b[Vv]alidation failed\b"),
    }

    # Sub-patterns for specific messages and TID extraction points
    sub_patterns = {
        "discarding": [
            (
                re.compile(r"under the minimum of (\d+)"),
                "Size under minimum",
                r'Discarding\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"internal stop codons found"),
                "Internal stop codons found",
                r'Invalid ORF\(s\) for\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"incorrect fusions of splice junctions"),
                "Incorrect fusions of splice junctions",
                r'Discarded\s+["\']?(\S+?)["\']?\s+because',
            ),
            (
                re.compile(r"Both UTR presents with a truncated ORF"),
                "Both UTR present with truncated ORF",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"defined UTRs but no CDS feature"),
                "Defined UTRs but no CDS feature",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"has a CDS .*straddles 2 different exons"),
                "CDS which straddles 2 different exons",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"cannot reverse the strand of a coding transcript"),
                "Cannot reverse strand of coding transcript",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"Overlapping exons found"),
                "Overlapping exons found",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"Invalid number of coding exons"),
                "Invalid number of coding exons",
                r'Invalid number of coding exons for\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"Invalid start and stop of the ORF"),
                "Invalid start and stop of the ORF",
                r'Discarded generically invalid transcript\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"debords its exon"),
                "Debords its exon",
                r'in\s+["\']?(\S+?)["\']?\s+debords',
            ),
            (
                re.compile(r"5\'UTR present with a truncated ORF"),
                "5'UTR present with a truncated ORF",
                r"Discarded generically invalid transcript ([^\s,]+)",
            ),
        ],
        "validation": [
            (
                re.compile(r"start must be less than end"),
                "Assertion failure start must be less than end",
                r'Validation failed on\s+["\']?(\S+?)["\']?[,\s]',
            ),
            (
                re.compile(r"assertion failure"),
                "Assertion failure",
                r'Validation failed on\s+["\']?(\S+?)["\']?[,\s]',
            ),
        ],
        "excluding": [
            (
                re.compile(r"\bExcluding\b.*\bas redundant with\b"),
                "Redundant",
                r'Excluding\s+["\']?(\S+?)["\']?\s+as redundant',
            )
        ],
    }


    # Function to classify each line
    def classify_line(line, stats, unique_tid_counts, log_entries):
        entry = {"Category": None, "SubCategory": None, "Message": line, "tid": None}
        tid = None

        # Check if line matches any main category pattern
        matched_category = None
        for category, pattern in main_patterns.items():
            if pattern.search(line):
                matched_category = category
                break

        # Skip processing if no main category pattern is matched
        if not matched_category:
            return

        # Process subcategory patterns if a category is matched
        for sub_pattern, description, tid_regex in sub_patterns[matched_category]:
            if sub_pattern.search(line):
                entry["Category"] = matched_category
                entry["SubCategory"] = description
                tid_match = re.search(tid_regex, line)
                if tid_match:
                    tid = tid_match.group(1).strip("\"' ,.")
                    entry["tid"] = tid
                key = (entry["Category"], entry["SubCategory"])
                stats[key] += 1
                if tid:
                    unique_tid_counts[key].add(tid)
                log_entries.append(entry)
                return

        # Categorize as "General" if no specific subcategory matched
        entry["Category"] = matched_category
        entry["SubCategory"] = "General"
        key = (entry["Category"], "General")
        stats[key] += 1  # Increment the count for General category
        log_entries.append(entry)


    # Data structure to collect plot data
    plot_data = defaultdict(lambda: defaultdict(int))

    # Process each log file
    for log_file_path in log_file_paths:
        if not os.path.isfile(log_file_path):
            print(f"Error: File '{log_file_path}' not found.")
            continue

        summary_output_dir = os.path.dirname(os.path.abspath(log_file_path))
        # Remove ".prepare.log" suffix if present from the file name for labeling
        base_name = os.path.splitext(os.path.basename(log_file_path))[0].replace(
            ".prepare", ""
        )

        stats = {key: 0 for key in all_combinations}
        unique_tid_counts = {key: set() for key in all_combinations}
        log_entries = []

        with open(log_file_path, "r") as file:
            for line in file:
                classify_line(line.strip(), stats, unique_tid_counts, log_entries)

        # Populate plot_data for plotting
        for key, value in stats.items():
            category, subcategory = key
            plot_data[base_name][subcategory] = value

        # Define output paths for parsed summary and summary statistics
        output_file_path = os.path.join(summary_output_dir, f"{base_name}_parsed_summary.csv")
        summary_file_path = os.path.join(summary_output_dir, f"{base_name}_summary_stats.csv")

        # Write parsed log entries to CSV
        with open(output_file_path, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["Category", "SubCategory", "tid", "Message"])
            for entry in log_entries:
                csv_writer.writerow(
                    [
                        entry["Category"],
                        entry["SubCategory"],
                        entry["tid"],
                        entry["Message"],
                    ]
                )

        # Calculate unique TID counts across all categories and subcategories
        all_unique_tids = set()
        for tids in unique_tid_counts.values():
            all_unique_tids.update(tids)

        # Write summary statistics to CSV
        with open(summary_file_path, "w", newline="") as summaryfile:
            summary_writer = csv.writer(summaryfile)
            summary_writer.writerow(
                ["Category", "SubCategory", "Count", "Unique TID Count"]
            )

            for key in all_combinations:
                count = stats.get(key, 0)
                unique_tid_count = len(unique_tid_counts.get(key, set()))
                summary_writer.writerow([key[0], key[1], count, unique_tid_count])

                # Collect data for plotting
                plot_data[base_name][key[1]] = count

            # Add totals for all categories
            total_count = sum(stats.values())
            total_unique_tid_count = len(all_unique_tids)
            summary_writer.writerow(
                ["Total", "All Categories", total_count, total_unique_tid_count]
            )

        print(f"Parsing complete for {log_file_path}.")
        print(f"Extracted log information saved to {output_file_path}")
        print(f"Summary stats saved to {summary_file_path}")

    # Convert plot_data to a DataFrame
    plot_df = pd.DataFrame(plot_data).fillna(0)

    # Generate stacked bar plot and save as PNG
    if plot_df.empty:
        print("No data to plot. Skipping plot generation.")
        sys.exit(0)

    labels = plot_df.columns.tolist()

    fig, ax = plt.subplots(figsize=(10, 9))

    color_palette = {
        "Size under minimum": "#1f77b4",
        "Internal stop codons found": "#ff7f0e",
        "Incorrect fusions of splice junctions": "#2ca02c",
        "Both UTR present with truncated ORF": "#d62728",
        "Defined UTRs but no CDS feature": "#9467bd",
        "CDS which straddles 2 different exons": "#8c564b",
        "Cannot reverse strand of coding transcript": "#e377c2",
        "Overlapping exons found": "#7f7f7f",
        "Invalid number of coding exons": "#bcbd22",
        "Invalid start and stop of the ORF": "#17becf",
        "Debords its exon": "#bcbddc",
        "5'UTR present with a truncated ORF": "#cc6666",
        "Assertion failure start must be less than end": "#ff9896",
        "Assertion failure": "#c5b0d5",
        "Redundant": "#c49c94",
        "General": "#5A5A5A",
    }

    bottoms = [0] * len(labels)

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

    plt.ylabel("Counts")
    plt.xlabel("Sample")
    plt.title(plot_title, loc="center", pad=20)
    plt.xticks(rotation=45, ha="right", fontsize=8)

    # Adjust legend to appear below the plot
    plt.legend(
        title="Categories",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=3,
        fontsize=8,
        title_fontsize=10,
    )
    plt.subplots_adjust(bottom=0.2)
    plt.tight_layout()

    # Save the plot using the specified output base name
    output_plot_file = f"{plot_output_base}_summary_plot.png"
    plt.savefig(output_plot_file)

    print(f"Plot saved to '{output_plot_file}'")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lifton/Liftoff transcript filtering and comparison pipeline

"""

import os
import sys
import argparse
import logging
import shutil

script = os.path.basename(sys.argv[0])
executed_command = " ".join(sys.argv)

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)


class HelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


# list of excluded categories from plot
EXCLUDED_CATEGORIES = [
    "\t5'UTR present with a truncated ORF",
    "Assertion failure",
    "Assertion failure start must be less than end",
    "Both UTR present with truncated ORF",
    "CDS which straddles 2 different exons",
    "Cannot reverse strand of coding transcript",
    "Debords its exon",
    "Defined UTRs but no CDS feature",
    "Duplicate parent feature ID",
    "General",
    "Incorrect fusions of splice junctions",
    "Internal stop codons found",
    "Invalid CDS length",
    "Invalid number of coding exons",
    "Invalid start and stop of the ORF",
    "Overlapping CDS",
    "Overlapping exons found",
    "Redundant",
    "Seqid mismatch*",
    "Short intron",
    "Size under minimum",
    "Strand conflict child*",
    "Strand conflict gene*",
    "Strand conflict gene-child*",
]

help_text_note = f"""
Note:
  - In --single mode, you must provide exactly one of --lifton_gff or --liftoff_gff.
  - In paired mode, both --lifton_gff and --liftoff_gff are required.
  - The script will symlink input files into the output directory.
  - Filtering is based on Mikado Prepare identified errors. Categories of errors can be excluded i.e. not filtered via --exclude_from_filtering. Full list of categories can be found below:
    {",\n\t".join(EXCLUDED_CATEGORIES)}
    * Categories marked with * cannot be excluded from being filtered.

   - The script assumes the presence of external tools:
      - gffread
      - gt (GenomeTools)
      - mikado
"""


class FilterLiftonLiftoff:
    def __init__(self, args):
        self.args = args
        self.output = str()
        self.analysis_dir = str()
        self.prefix = None
        self.force = self.args.force
        self.debug = self.args.debug
        self.lifton_label = None
        self.liftoff_label = None
        self.did_lifton = False
        self.did_liftoff = False
        self.lifton_sorted_gff = None
        self.liftoff_sorted_gff = None
        self.genome = None
        self.lifton_gff = None
        self.liftoff_gff = None
        self.mikado_list_label = "LiftClean_PREFIX"

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)
            logging.debug("Debug mode enabled")

        # check if input and output files/directories exist
        if not os.path.isfile(self.args.genome_fasta):
            logging.error(f"Reference genome file {self.args.genome_fasta} not found")
            sys.exit(1)
        if self.args.lifton_gff and not os.path.isfile(self.args.lifton_gff):
            logging.error(f"Lifton GFF file {self.args.lifton_gff} not found")
            sys.exit(1)
        if self.args.liftoff_gff and not os.path.isfile(self.args.liftoff_gff):
            logging.error(f"Liftoff GFF file {self.args.liftoff_gff} not found")
            sys.exit(1)
        # remove output directory if it exists and --force is enabled
        if os.path.isdir(self.args.output) and self.force:
            logging.info(
                f"Output directory {self.args.output} exists, removing it due to --force option"
            )
            shutil.rmtree(self.args.output)

        if self.args.output and not os.path.isdir(self.args.output):
            logging.info(f"Output directory {self.args.output} not found, creating it")
            os.makedirs(self.args.output, exist_ok=True)

        self.output = os.path.abspath(self.args.output)

        self.analysis_dir = os.path.join(self.output, "analysis")
        if not os.path.isdir(self.analysis_dir):
            logging.info(f"Creating analysis directory {self.analysis_dir}")
            os.makedirs(self.analysis_dir, exist_ok=True)

        if self.args.single:
            logging.info("Single GFF mode enabled")
            # only one GFF file allowed, exit if both or none provided
            if self.args.lifton_gff and self.args.liftoff_gff:
                logging.error(
                    "Both Lifton and Liftoff GFF files provided, but --single option enabled. Please provide only one GFF file."
                )
            # When --single is enabled, only one GFF file is allowed; exit if both or none are provided
            if not self.args.lifton_gff and not self.args.liftoff_gff:
                logging.error(
                    "No GFF file provided, but --single option enabled. Please provide one GFF file."
                )
                sys.exit(1)

        # enable self.args.single mode even if --single is not provided but only one GFF is given
        if not self.args.single:
            if self.args.lifton_gff and not self.args.liftoff_gff:
                logging.warning("Only Lifton GFF file provided, enabling --single mode")
                self.args.single = True
            elif self.args.liftoff_gff and not self.args.lifton_gff:
                logging.warning(
                    "Only Liftoff GFF file provided, enabling --single mode"
                )
                self.args.single = True
            elif not self.args.lifton_gff and not self.args.liftoff_gff:
                logging.error("No GFF files provided. Please provide at least one.")
                sys.exit(1)

        if self.args.alt_lifton_label:
            self.lifton_label = self.args.alt_lifton_label
        if self.args.alt_liftoff_label:
            self.liftoff_label = self.args.alt_liftoff_label
        if self.args.prefix:
            if self.lifton_label:
                self.lifton_label = f"{self.lifton_label}_{self.args.prefix}"
            else:
                self.lifton_label = f"lifton_{self.args.prefix}"
            if self.liftoff_label:
                self.liftoff_label = f"{self.liftoff_label}_{self.args.prefix}"
            else:
                self.liftoff_label = f"liftoff_{self.args.prefix}"

        # link the input files to output directory
        # - genome as genome.fasta
        # - lifton_gff as lifton.gff
        # - liftoff_gff as liftoff.gff
        genome_link = os.path.join(self.output, "genome.fasta")
        if not os.path.isfile(genome_link):
            org_file = os.path.abspath(self.args.genome_fasta)
            os.symlink(org_file, genome_link)
            logging.info(f"Linked genome file {org_file} to {genome_link}")
        self.genome = os.path.join(self.output, "genome.fasta")
        if self.args.lifton_gff:
            lifton_link = os.path.join(self.output, f"{self.lifton_label}.gff")
            if not os.path.isfile(lifton_link):
                org_file = os.path.abspath(self.args.lifton_gff)
                os.symlink(org_file, lifton_link)
                logging.info(f"Linked Lifton GFF file {org_file} to {lifton_link}")
            self.lifton_gff = os.path.join(self.output, f"{self.lifton_label}.gff")
        if self.args.liftoff_gff:
            liftoff_link = os.path.join(self.output, f"{self.liftoff_label}.gff")
            if not os.path.isfile(liftoff_link):
                org_file = os.path.abspath(self.args.liftoff_gff)
                os.symlink(org_file, liftoff_link)
                logging.info(f"Linked Liftoff GFF file {org_file} to {liftoff_link}")
            self.liftoff_gff = os.path.join(self.output, f"{self.liftoff_label}.gff")

        if self.args.prefix:
            self.prefix = self.args.prefix
        else:
            self.prefix = os.path.basename(os.path.abspath(self.args.output))
        logging.info(f"Using prefix '{self.prefix}' for plot titles")
        logging.info(f"Output will be saved to {self.output}")

    def process_gff(self, gff_file, label):
        # strand checks
        corrected = os.path.join(self.analysis_dir, f"{label}.corrected.gff")
        corrected_log = os.path.join(
            self.analysis_dir, f"{label}.corrected.gff_strand_checker.log"
        )

        # intron checks
        intron_corrected = os.path.join(
            self.analysis_dir, f"{label}.short_introns.corrected.gff"
        )
        intron_corrected_tsv = os.path.join(
            self.analysis_dir, f"{label}.short_introns.tsv"
        )
        intron_corrected_mapping_tsv = os.path.join(
            self.analysis_dir, f"{label}.short_introns.mapping.tsv"
        )

        # gffread
        corrected_gffread = os.path.join(
            self.analysis_dir, f"{label}.corrected.gffread.gff"
        )

        # gt gff3
        sorted_gff = os.path.join(self.analysis_dir, f"{label}.sorted.gff")
        sorted_gff_log = os.path.join(self.analysis_dir, f"{label}.sorted.gff.log")

        # mikado prepare
        mikado_dir = os.path.join(self.analysis_dir, f"mikado_prepare_{label}")
        mikado_log = os.path.join(mikado_dir, f"mikado_prepare_{label}.log")

        # strand checks
        logging.info(f"Strand-checking {label} GFF file")
        cmd = "gff_strand_checker --check_parent_seq_id "
        if self.args.check_dup_ids:
            cmd += "--check_dup_ids "
        if self.args.allow_trans_splicing:
            cmd += "--allow_trans_splicing "
        cmd += (
            f" --remove {gff_file} --output_dir {self.analysis_dir} > {corrected_log}"
        )

        if not os.path.isfile(corrected):
            self.process_cmd(cmd)
        if os.path.isfile(corrected) and self.debug:
            logging.debug(f"Skipping strand-checking, {corrected} already exists")
        if not os.path.isfile(corrected):
            logging.error(f"Corrected GFF file {corrected} not created")
            sys.exit(1)

        # intron checks
        logging.info(
            f"Filtering models with introns <= {self.args.min_intron_length} bp"
        )
        cmd = (
            "filter_short_introns "
            f"--intron_size {self.args.min_intron_length} "
            f"--intron_size_ends {self.args.min_intron_length} "
            f"--detailed_discard_list_file {intron_corrected_tsv} "
            f"--discard_list_mapping_file {intron_corrected_mapping_tsv} "
        )
        if self.args.transcript_types:
            cmd += f"--transcript_types {self.args.transcript_types} "
        if self.args.gene_types:
            cmd += f"--gene_types {self.args.gene_types} "
        cmd += f"{corrected} --output {intron_corrected}"

        if not os.path.isfile(intron_corrected):
            self.process_cmd(cmd)
        if os.path.isfile(intron_corrected) and self.debug:
            logging.debug(
                f"Skipping intron filtering, {intron_corrected} already exists"
            )
        if not os.path.isfile(intron_corrected):
            logging.error(f"Intron-corrected GFF file {intron_corrected} not created")
            sys.exit(1)

        # gffread
        logging.info(f"Running gffread on {label} GFF file")
        cmd = f"gffread {self.args.gffread_params} -o {corrected_gffread} {intron_corrected}"
        if not os.path.isfile(corrected_gffread):
            self.process_cmd(cmd)
        if os.path.isfile(corrected_gffread) and self.debug:
            logging.debug(f"Skipping gffread, {corrected_gffread} already exists")
        if not os.path.isfile(corrected_gffread):
            logging.error(f"Corrected GFFREAD file {corrected_gffread} not created")
            sys.exit(1)

        # gt gff3
        logging.info(f"Sorting {label} GFF file")
        cmd = f"gt gff3 {self.args.gt_gff3_params} {corrected_gffread} > {sorted_gff} 2> {sorted_gff_log}"
        if not os.path.isfile(sorted_gff):
            self.process_cmd(cmd)
        if os.path.isfile(sorted_gff) and self.debug:
            logging.debug(f"Skipping sorting, {sorted_gff} already exists")
        if not os.path.isfile(sorted_gff) or os.path.getsize(sorted_gff) == 0:
            logging.error(f"Sorted GFF file {sorted_gff} not created or empty")
            sys.exit(1)

        # mikado prepare
        logging.info(
            f"Running mikado prepare on {label} GFF file (min cDNA length {self.args.minimum_cdna_length})"
        )

        # create list file
        if not os.path.isdir(mikado_dir):
            os.makedirs(mikado_dir, exist_ok=True)
        list_file = os.path.join(mikado_dir, "list.txt")
        with open(list_file, "w") as f:
            if self.args.limit_filters:
                # marking the model as reference 'True' to avoid filtering by mikado prepare
                f.write(f"{sorted_gff}\t{self.mikado_list_label}\tTrue\t0\tTrue\n")
            else:
                f.write(f"{sorted_gff}\t{self.mikado_list_label}\tTrue\t0\tFalse\n")
        logging.info(f"Created Mikado list file {list_file}")

        cmd = (
            "mikado prepare "
            f"--fasta {self.genome} "
            f"--minimum-cdna-length {self.args.minimum_cdna_length} "
            f"-p {self.args.threads} "
            f"-od {mikado_dir} "
            f"-o mikado_prepare_{label}.gtf "
            f"-of mikado_prepare_{label}.fasta "
            f"-l {mikado_log} "
            f"{self.args.mikado_prepare_params}"
            f"--list {list_file}"
        )
        if not os.path.isfile(os.path.join(mikado_dir, f"mikado_prepare_{label}.gtf")):
            self.process_cmd(cmd)
        if (
            os.path.isfile(os.path.join(mikado_dir, f"mikado_prepare_{label}.gtf"))
            and self.debug
        ):
            logging.debug(
                f"Skipping mikado prepare, {os.path.join(mikado_dir, f'mikado_prepare_{label}.gtf')} already exists"
            )
        if not os.path.isfile(os.path.join(mikado_dir, f"mikado_prepare_{label}.gtf")):
            logging.error(
                f"Mikado GTF file {os.path.join(mikado_dir, f'mikado_prepare_{label}.gtf')} not created"
            )
            sys.exit(1)

        logging.info(f"Finished processing {label} GFF file")

        if label == self.lifton_label:
            self.lifton_sorted_gff = sorted_gff
            self.did_lifton = True
        elif label == self.liftoff_label:
            self.liftoff_sorted_gff = sorted_gff
            self.did_liftoff = True

    def generate_summary_stats(self):
        """
        Generate summary_stats.csv via parse_prepare_log
        """
        rm_prefix = self.mikado_list_label + "_"
        if self.args.single:
            if self.did_lifton:
                cmd = (
                    "parse_prepare_log "
                    f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}.log')} "
                    f"--title 'Rejected Transcripts {self.lifton_label}' "
                    f"--rm_prefix {rm_prefix} "
                    f"--output_prefix '{self.prefix}_Rejected_Transcripts' "
                    f"--output {self.output}"
                )
                self.process_cmd(cmd)
            else:
                cmd = (
                    "parse_prepare_log "
                    f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}.log')} "
                    f"--title 'Rejected Transcripts {self.liftoff_label}' "
                    f"--rm_prefix {rm_prefix} "
                    f"--output_prefix '{self.prefix}_Rejected_Transcripts' "
                    f"--output {self.output}"
                )
                self.process_cmd(cmd)
        else:
            cmd = (
                "parse_prepare_log "
                f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}.log')} "
                f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}.log')} "
                f"--title 'Rejected Transcripts {self.lifton_label} and {self.liftoff_label}' "
                f"--rm_prefix {rm_prefix} "
                f"--output_prefix '{self.prefix}_Rejected_Transcripts' "
                f"--output {self.output}"
            )
            self.process_cmd(cmd)

    def append_computed_summary(self, label):
        """
        Append strand-check and intron filtering summary
        """
        logging.info(
            f"Appending strand-check and intron filtering summary for {label} to mikado_prepare_{label}_summary_stats.csv"
        )

        summary_stats_csv = os.path.join(
            self.analysis_dir,
            f"mikado_prepare_{label}",
            f"mikado_prepare_{label}_summary_stats.csv",
        )
        strand_corrected_log = os.path.join(
            self.analysis_dir, f"{label}.corrected.gff_strand_checker.log"
        )
        intron_corrected_tsv = os.path.join(
            self.analysis_dir, f"{label}.short_introns.tsv"
        )

        strand_append = False
        new_summary_lines = []
        # check whether strand_corrected_log has '^Category,SubCategory'
        # if yes, skip header and append all lines until '^Checks complete' to list
        with open(strand_corrected_log, "r") as f:
            for line in f:
                if line.startswith("Category,SubCategory"):
                    strand_append = True
                    continue
                if line.startswith("Checks complete"):
                    strand_append = False
                    break
                if strand_append:
                    new_summary_lines.append(line)

        intron_append = False
        # check whether intron_corrected_tsv has '^Category,SubCategory'
        # if yes, skip header and append all lines until '^Checks complete' to list
        with open(intron_corrected_tsv, "r") as f:
            for line in f:
                if line.startswith("Category,SubCategory"):
                    intron_append = True
                    continue
                if line.startswith("Checks complete"):
                    intron_append = False
                    break
                if intron_append:
                    new_summary_lines.append(line)

        # Save the original CSV body, skip header and drop existing Total, then recalculate totals
        if new_summary_lines:
            logging.info(
                f"Prepending strand-check and intron filter summary for {label} from {strand_corrected_log} and {intron_corrected_tsv} to {summary_stats_csv}"
            )

            with open(summary_stats_csv, "r") as f:
                original_lines = f.readlines()
            header = original_lines[0]
            body_lines = [
                line for line in original_lines[1:] if not line.startswith("Total,")
            ]
            merged_lines = [header] + body_lines + new_summary_lines
            # Recalculate totals (sum columns 3 & 4)
            total_count = sum(
                int(line.split(",")[2]) for line in merged_lines[1:] if line.strip()
            )
            total_uid = sum(
                int(line.split(",")[3]) for line in merged_lines[1:] if line.strip()
            )
            merged_lines.append(f"Total,All Categories,{total_count},{total_uid}\n")
            with open(summary_stats_csv, "w") as f:
                f.writelines(merged_lines)
                f.write("\n")
        else:
            logging.warning(
                f"No strand-check and intron filter summary lines found in {strand_corrected_log} and {intron_corrected_tsv}, skipping append"
            )

    def regenerate_summary_plot(self):
        """
        Regenerating summary plot from updated summary_stats.csv
        """
        if self.args.single:
            if self.did_lifton:
                cmd = (
                    "plot_summary_csvs "
                    f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}_summary_stats.csv')} "
                    f"--title 'Rejected Transcripts {self.lifton_label}' "
                    f"--output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' "
                    f"--output {self.output}"
                )
                self.process_cmd(cmd)
            else:
                cmd = (
                    "plot_summary_csvs "
                    f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}_summary_stats.csv')} "
                    f"--title 'Rejected Transcripts {self.liftoff_label}' "
                    f"--output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' "
                    f"--output {self.output}"
                )
                self.process_cmd(cmd)
        else:
            cmd = (
                "plot_summary_csvs "
                f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}_summary_stats.csv')} "
                f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}_summary_stats.csv')} "
                f"--title 'Comparison of Rejected Transcripts {self.prefix}' "
                f"--output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' "
                f"--output {self.output}"
            )
            self.process_cmd(cmd)

    def compare_parsed_summary_csv(self):
        """
        Running compare_parsed_summary_csv.py (single/paired mode)
        """
        if self.args.single:
            logging.info("Running compare_parsed_summary_csv.py (single)")
            if self.did_lifton:
                cmd = (
                    "compare_parsed_summary_csv "
                    "--single "
                    f"--gff_file {self.lifton_sorted_gff} "
                    f"--csv_file {os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}_parsed_summary.csv')} "
                    f"--output {self.output} "
                    f"--exclude '{self.args.exclude_from_filtering}'"
                )
                self.process_cmd(cmd)
            else:
                cmd = (
                    "compare_parsed_summary_csv "
                    "--single "
                    f"--gff_file {self.liftoff_sorted_gff} "
                    f"--csv_file {os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}_parsed_summary.csv')} "
                    f"--output {self.output} "
                    f"--exclude '{self.args.exclude_from_filtering}'"
                )
                self.process_cmd(cmd)
        else:
            logging.info("Running compare_parsed_summary_csv.py (paired)")
            cmd = (
                "compare_parsed_summary_csv "
                f"--gff_files {self.lifton_sorted_gff} {self.liftoff_sorted_gff} "
                f"--csv_files {os.path.join(self.analysis_dir, f'mikado_prepare_{self.lifton_label}', f'mikado_prepare_{self.lifton_label}_parsed_summary.csv')} "
                f"{os.path.join(self.analysis_dir, f'mikado_prepare_{self.liftoff_label}', f'mikado_prepare_{self.liftoff_label}_parsed_summary.csv')} "
                f"--labels {self.lifton_label},{self.liftoff_label} "
                f"--output_prefix {self.prefix} "
                f"--output {self.output} "
                f"--plot_title 'Comparison of Retained and Rejected IDs {self.prefix}' "
                f"--exclude '{self.args.exclude_from_filtering}'"
            )
            self.process_cmd(cmd)

    def extract_gene_id_gff(self, label):
        """
        Running extract_gene_id_gff for retained IDs
        """
        id_file = os.path.join(self.analysis_dir, f"{label}.sorted.ids.retained.txt")
        out_tsv = os.path.join(self.output, f"{label}.sorted.ids.retained.tsv")
        logging.info(f"Running extract_gene_id_gff for {label}")
        cmd = None
        if label == self.lifton_label:
            cmd = f"extract_gene_id_gff {id_file} {self.lifton_sorted_gff} {out_tsv}"
        else:
            cmd = f"extract_gene_id_gff {id_file} {self.liftoff_sorted_gff} {out_tsv}"
        self.process_cmd(cmd)

    def mikado_util_grep(self, label):
        """
        Running mikado util grep
        """
        tsv_file = os.path.join(self.output, f"{label}.sorted.ids.retained.tsv")
        out_gff = os.path.join(self.output, f"{label}.sorted.retained.gff")
        logging.info(f"Running mikado util grep for {label}")
        cmd = None
        if label == self.lifton_label:
            cmd = f"mikado util grep {tsv_file} {self.lifton_sorted_gff} {out_gff}"
        else:
            cmd = f"mikado util grep {tsv_file} {self.liftoff_sorted_gff} {out_gff}"
        self.process_cmd(cmd)

    def final_summary(self):
        """
        Final summary of retained files
        """
        print("\n=== ✅ Pipeline Completed Successfully ===")
        if self.did_lifton:
            print(
                f" - lifton retained GFF: {os.path.join(self.output, f'{self.lifton_label}.sorted.retained.gff')}"
            )
            print(
                f" - lifton retained TSV: {os.path.join(self.output, f'{self.lifton_label}.sorted.ids.retained.tsv')}"
            )
        if self.did_liftoff:
            print(
                f" - liftoff retained GFF: {os.path.join(self.output, f'{self.liftoff_label}.sorted.retained.gff')}"
            )
            print(
                f" - liftoff retained TSV: {os.path.join(self.output, f'{self.liftoff_label}.sorted.ids.retained.tsv')}"
            )
        comparison_csv = os.path.join(
            self.output, f"{self.prefix}_Rejected_Transcripts.csv"
        )
        if os.path.isfile(comparison_csv):
            print(f" - Comparison CSV: {comparison_csv}")
        print()

    def process_cmd(self, cmd):
        # capture exit code and log error if non-zero
        logging.info(f"Running command: {cmd}")
        exit_code = os.system(cmd)
        if exit_code != 0:
            logging.error(f"Command failed with exit code {exit_code}")
            sys.exit(exit_code)

    def run(self):
        logging.info("Starting LiftClean")
        logging.info(f"Command: {executed_command}")

        #######################################
        #    process_gff:
        #    strand-check → gffread → sort → mikado prepare
        #######################################
        if self.lifton_gff:
            logging.info("Processing Lifton GFF file")
            self.process_gff(self.lifton_gff, label=self.lifton_label)
        if self.liftoff_gff:
            logging.info("Processing Liftoff GFF file")
            self.process_gff(self.liftoff_gff, label=self.liftoff_label)

        #######################################
        #    Generate summary_stats.csv via parse_prepare_log
        #######################################
        self.generate_summary_stats()

        #######################################
        #    Prepend strand-check summary & recalc Total
        #######################################
        if self.did_lifton:
            self.append_computed_summary(label=self.lifton_label)
        if self.did_liftoff:
            self.append_computed_summary(label=self.liftoff_label)

        #######################################
        #    Regenerate summary plot
        #######################################
        self.regenerate_summary_plot()

        #######################################
        #    compare_parsed_summary_csv.py
        #######################################
        self.compare_parsed_summary_csv()

        #######################################
        #    extract_gene_id_gff.py
        #######################################
        if self.lifton_gff:
            self.extract_gene_id_gff(label=self.lifton_label)
        if self.liftoff_gff:
            self.extract_gene_id_gff(label=self.liftoff_label)

        #######################################
        #    mikado util grep
        #######################################
        if self.lifton_gff:
            self.mikado_util_grep(label=self.lifton_label)
        if self.liftoff_gff:
            self.mikado_util_grep(label=self.liftoff_label)

        #######################################
        #    Final summary
        #######################################
        self.final_summary()

        # Here would be the main logic of the filtering and comparison
        logging.info("Finished LiftClean")


def main():
    parser = argparse.ArgumentParser(
        prog=script,
        formatter_class=HelpFormatter,
        description="""
        Lifton/Liftoff transcript filtering and comparison pipeline
        """,
        epilog=f"{help_text_note}\n",
    )
    parser.add_argument(
        "-g",
        "--genome_fasta",
        required=True,
        help="Provide reference genome FASTA file",
    )
    # --lifton_gff
    parser.add_argument(
        "-n",
        "--lifton_gff",
        help="Provide Lifton output GFF file",
    )
    # --liftoff_gff
    parser.add_argument(
        "-f",
        "--liftoff_gff",
        help="Provide Liftoff output GFF file",
    )
    # alternative lifton label
    parser.add_argument(
        "--alt_lifton_label",
        type=str,
        default="lifton",
        help="Alternative label for Lifton in plots and outputs [default:%(default)s]",
    )
    # alternative liftoff label
    parser.add_argument(
        "--alt_liftoff_label",
        type=str,
        default="liftoff",
        help="Alternative label for Liftoff in plots and outputs [default:%(default)s]",
    )
    # --prefix for plots titles
    parser.add_argument(
        "-p",
        "--prefix",
        help="Provide a label for file name prefix and plot titles to distinguish images from multiple runs. If provided, it will be suffixed to the --alt_lifton_label and --alt_liftoff_label. For example, you can use 'Wheat_Accession1' and the output will be labeled accordingly like lifton_Wheat_Accession1 or liftoff_Wheat_Accession1 [default:%(default)s]",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="Provide output directory [default:%(default)s]",
    )
    parser.add_argument(
        "-s",
        "--single",
        action="store_true",
        help="Only one GFF allowed. Default is to process both --lifton_gff and --liftoff_gff inputs [default:%(default)s]",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        help="Number of threads to use for copying",
    )
    parser.add_argument(
        "-e",
        "--exclude_from_filtering",
        type=str,
        default="",
        help="Provide a comma-separated list of types to be excluded from being filtered (default: none, e.g. 'Size under minimum,Redundant' etc [default:%(default)s]",
    )
    # --minimum-cdna-length
    parser.add_argument(
        "-m",
        "--minimum_cdna_length",
        type=int,
        default=48,
        help="Provide minimum cDNA length for filtering. Anything lower will be excluded [default: %(default)s]",
    )
    parser.add_argument(
        "-i",
        "--min_intron_length",
        type=int,
        default=0,
        help="Provide intron size. Any models with intron size under this value will be removed. Bookended and overlapping exon features are merged with the default setting [default: %(default)s]",
    )
    parser.add_argument(
        "--limit_filters",
        action="store_true",
        help="Enable this to not exclude models based on categories: Size under minimum, Incorrect fusions of splice junctions, Cannot reverse strand of coding transcript, Redundant [default:%(default)s]",
    )
    parser.add_argument(
        "--check_dup_ids",
        action="store_false",
        help="Check for duplicate IDs, but only for features that are parents; if found, remove all instances and all descendants [default:%(default)s]",
    )
    parser.add_argument(
        "--allow_trans_splicing",
        action="store_true",
        help="If set, skip all checks for any parent (and its children) where the parent has 'exception=trans-splicing' [default:%(default)s]",
    )
    # Configurable transcript and gene types
    parser.add_argument(
        "--transcript_types",
        default="mRNA,primary_transcript,transcript,lnc_RNA,ncRNA,miRNA,rRNA,tRNA,snoRNA,snRNA,scaRNA,pseudogenic_transcript,antisense_RNA",
        help="Comma-separated list of feature types to treat as transcripts (default: %(default)s)",
    )
    parser.add_argument(
        "--gene_types",
        default="gene,ncRNA_gene,pseudogene",
        help="Comma-separated list of feature types to treat as gene-level (default: %(default)s)",
    )
    # add gffread params
    parser.add_argument(
        "--gffread_params",
        default="--keep-genes -F",
        type=str,
        help='Additional parameters to pass to gffread (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --gffread_params="<params>" [default:"%(default)s"]',
    )
    # gt gff3 params
    parser.add_argument(
        "--gt_gff3_params",
        default="-sort -tidy -retainids yes",
        type=str,
        help='Additional parameters to pass to gt (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --gt_gff3_params="<params>" [default:"%(default)s"]',
    )
    # mikado prepare params
    parser.add_argument(
        "--mikado_prepare_params",
        default="",
        type=str,
        help='Additional parameters to pass to mikado (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --mikado_prepare_params="<params>" [default:"%(default)s"]',
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Enable this option to force removal of pre-existing output folder [default:%(default)s]",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable this option for debugging [default:%(default)s]",
    )
    args = parser.parse_args()

    # check if gffread, gt, mikado are available
    logging.info("Checking for required tools in PATH")
    for tool in [
        "gffread",
        "gt",
        "mikado",
    ]:
        if not shutil.which(tool):
            logging.error(f"Required tool '{tool}' not found in PATH")
            sys.exit(1)
        else:
            logging.info(
                f"Found required tool in PATH: '{tool}' - {shutil.which(tool)}"
            )

    FilterLiftonLiftoff(args).run()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE

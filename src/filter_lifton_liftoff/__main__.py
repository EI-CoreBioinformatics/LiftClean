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

__author__ = "David Swarbreck"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "David.Swarbreck@earlham.ac.uk"

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
    "\tSize under minimum",
    "Internal stop codons found",
    "Incorrect fusions of splice junctions",
    "Both UTR present with truncated ORF",
    "Defined UTRs but no CDS feature",
    "CDS which straddles 2 different exons",
    "Cannot reverse strand of coding transcript",
    "Overlapping exons found",
    "Invalid number of coding exons",
    "Invalid start and stop of the ORF",
    "Debords its exon",
    "5'UTR present with a truncated ORF",
    "Assertion failure start must be less than end",
    "Assertion failure",
    "Redundant",
    "Seqid mismatch",
    "Strand conflict child",
    "Strand conflict gene",
]

help_text_note = f"""
Note:
  - In --single mode, you must provide exactly one of --lifton_gff or --liftoff_gff.
  - In paired mode, both --lifton_gff and --liftoff_gff are required.
  - The script will symlink input files into the output directory, the output directory is also used as a prefix.
  - Filtering is based on Mikado Prepare identified errors, categories of errors can be excluded i.e. not filtered via --exclude. Full list of categories can be found below:
    {',\n\t'.join(EXCLUDED_CATEGORIES)}
  - The script assumes the presence of external tools:
      - gffread
      - gt (GenomeTools)
      - mikado
"""


class FilterLiftonLiftoff:
    def __init__(self, args):
        self.args = args
        self.output = None
        self.analysis_dir = None
        self.prefix = None
        self.force = self.args.force
        self.debug = self.args.debug
        self.did_lifton = False
        self.did_liftoff = False
        self.lifton_sorted_gff = None
        self.liftoff_sorted_gff = None
        self.genome = None
        self.lifton_gff = None
        self.liftoff_gff = None

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
            lifton_link = os.path.join(self.output, "lifton.gff")
            if not os.path.isfile(lifton_link):
                org_file = os.path.abspath(self.args.lifton_gff)
                os.symlink(org_file, lifton_link)
                logging.info(f"Linked Lifton GFF file {org_file} to {lifton_link}")
            self.lifton_gff = os.path.join(self.output, "lifton.gff")
        if self.args.liftoff_gff:
            liftoff_link = os.path.join(self.output, "liftoff.gff")
            if not os.path.isfile(liftoff_link):
                org_file = os.path.abspath(self.args.liftoff_gff)
                os.symlink(org_file, liftoff_link)
                logging.info(f"Linked Liftoff GFF file {org_file} to {liftoff_link}")
            self.liftoff_gff = os.path.join(self.output, "liftoff.gff")

        if self.args.prefix:
            self.prefix = self.args.prefix
        else:
            self.prefix = os.path.basename(os.path.abspath(self.args.output))
        logging.info(f"Using prefix '{self.prefix}' for plot titles")
        logging.info(f"Output will be saved to {self.output}")

    def process_gff(self, gff_file, label):
        corrected = os.path.join(self.analysis_dir, f"{label}.corrected.gff")
        corrected_log = os.path.join(self.analysis_dir, f"{label}.corrected.gff_strand_checker.log")
        corrected_gffread = os.path.join(self.analysis_dir, f"{label}.corrected.gffread.gff")
        sorted_gff = os.path.join(self.analysis_dir, f"{label}.sorted.gff")
        sorted_gff_log = os.path.join(self.analysis_dir, f"{label}.sorted.gff.log")
        mikado_dir = os.path.join(self.analysis_dir, f"mikado_prepare_{label}")
        mikado_log = os.path.join(mikado_dir, f"mikado_prepare_{label}.log")

        logging.info(f"Strand-checking {label} GFF file")
        cmd = f"gff_strand_checker --check_parent_seq_id --remove {gff_file} --output_dir {self.analysis_dir} > {corrected_log}"
        if not os.path.isfile(corrected) or self.debug:
            self.process_cmd(cmd)
        else:
            logging.info(f"Skipping strand-checking, {corrected} already exists")
        if not os.path.isfile(corrected):
            logging.error(f"Corrected GFF file {corrected} not created")
            sys.exit(1)

        logging.info(f"Running gffread on {label} GFF file")
        cmd = f"gffread --keep-genes -Z -F -o {corrected_gffread} {corrected}"
        if not os.path.isfile(corrected_gffread) or self.debug:
            self.process_cmd(cmd)
        else:
            logging.info(f"Skipping gffread, {corrected_gffread} already exists")
        if not os.path.isfile(corrected_gffread):
            logging.error(f"Corrected GFFREAD file {corrected_gffread} not created")
            sys.exit(1)

        logging.info(f"Sorting {label} GFF file")
        cmd = f"gt gff3 -sort -tidy -retainids yes {corrected_gffread} > {sorted_gff} 2> {sorted_gff_log}"
        if not os.path.isfile(sorted_gff) or self.debug:
            self.process_cmd(cmd)
        else:
            logging.info(f"Skipping sorting, {sorted_gff} already exists")
        if not os.path.isfile(sorted_gff) or os.path.getsize(sorted_gff) == 0:
            logging.error(f"Sorted GFF file {sorted_gff} not created or empty")
            sys.exit(1)

        logging.info(
            f"Running mikado prepare on {label} GFF file (min cDNA length {self.args.minimum_cdna_length})"
        )
        cmd = f"mikado prepare --fasta {self.genome} --minimum-cdna-length {self.args.minimum_cdna_length} -p {self.args.threads} -od {mikado_dir} -o {f'mikado_prepare_{label}.gtf'} -of {f'mikado_prepare_{label}.fasta'} -l {mikado_log} {sorted_gff}"
        if (
            not os.path.isfile(os.path.join(mikado_dir, f"mikado_prepare_{label}.gtf"))
            or self.debug
        ):
            self.process_cmd(cmd)
        else:
            logging.info(
                f"Skipping mikado prepare, {os.path.join(mikado_dir, f'mikado_prepare_{label}.gtf')} already exists"
            )
        if not os.path.isfile(os.path.join(mikado_dir, f"mikado_prepare_{label}.gtf")):
            logging.error(
                f"Mikado GTF file {os.path.join(mikado_dir, f'mikado_prepare_{label}.gtf')} not created"
            )
            sys.exit(1)

        logging.info(f"Finished processing {label} GFF file")

        if label == "lifton":
            self.lifton_sorted_gff = sorted_gff
            self.did_lifton = True
        elif label == "liftoff":
            self.liftoff_sorted_gff = sorted_gff
            self.did_liftoff = True

    def generate_summary_stats(self):
        """
        Generate summary_stats.csv via parse_prepare_log
        """
        if self.args.single:
            if self.did_lifton:
                cmd = f"parse_prepare_log {os.path.join(self.analysis_dir, 'mikado_prepare_lifton', 'mikado_prepare_lifton.log')} --title 'Rejected Transcripts lifton' --output_prefix '{self.prefix}_Rejected_Transcripts' --output {self.output}"
                self.process_cmd(cmd)
            else:
                cmd = f"parse_prepare_log {os.path.join(self.analysis_dir, 'mikado_prepare_liftoff', 'mikado_prepare_liftoff.log')} --title 'Rejected Transcripts liftoff' --output_prefix '{self.prefix}_Rejected_Transcripts' --output {self.output}"
                self.process_cmd(cmd)
        else:
            cmd = f"parse_prepare_log {os.path.join(self.analysis_dir, 'mikado_prepare_lifton', 'mikado_prepare_lifton.log')} {os.path.join(self.analysis_dir, 'mikado_prepare_liftoff', 'mikado_prepare_liftoff.log')} --title 'Rejected Transcripts lifton and liftoff' --output_prefix '{self.prefix}_Rejected_Transcripts' --output {self.output}"
            self.process_cmd(cmd)

    def append_strand_summary(self, label):
        """
        Append strand-check summary helper
        """
        logging.info(
            f"Appending strand-check summary for {label} to mikado_prepare_{label}_summary_stats.csv"
        )

        summary_stats_csv = os.path.join(
            self.analysis_dir,
            f"mikado_prepare_{label}",
            f"mikado_prepare_{label}_summary_stats.csv",
        )
        corrected_log = os.path.join(self.analysis_dir, f"{label}.corrected.gff_strand_checker.log")

        start_append = False
        new_summary_lines = []
        # check whether corrected_log has '^Category,SubCategory'
        # if yes, skip header and append all lines until '^Checks complete' to list
        with open(corrected_log, "r") as f:
            for line in f:
                if line.startswith("Category,SubCategory"):
                    start_append = True
                    continue
                if line.startswith("Checks complete"):
                    start_append = False
                    break
                if start_append:
                    new_summary_lines.append(line)

        # Save the original CSV body, skip header and drop existing Total, then recalculate totals
        if new_summary_lines:
            logging.info(
                f"Prepending strand-check summary for {label} from {corrected_log} to {summary_stats_csv}"
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
                f"No strand-check summary lines found in {corrected_log}, skipping append"
            )

    def regenerate_summary_plot(self):
        """
        Regenerating summary plot from updated summary_stats.csv
        """
        if self.args.single:
            if self.did_lifton:
                cmd = f"plot_summary_csvs { os.path.join(self.analysis_dir, f"mikado_prepare_lifton", f"mikado_prepare_lifton_summary_stats.csv")} --title 'Rejected Transcripts lifton {self.prefix}' --output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' --output {self.output}"
                self.process_cmd(cmd)
            else:
                cmd = f"plot_summary_csvs { os.path.join(self.analysis_dir, f"mikado_prepare_liftoff", f"mikado_prepare_liftoff_summary_stats.csv")} --title 'Rejected Transcripts liftoff {self.prefix}' --output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' --output {self.output}"
                self.process_cmd(cmd)
        else:
            cmd = f"plot_summary_csvs {os.path.join(self.analysis_dir, 'mikado_prepare_lifton', 'mikado_prepare_lifton_summary_stats.csv')} {os.path.join(self.analysis_dir, 'mikado_prepare_liftoff', 'mikado_prepare_liftoff_summary_stats.csv')} --title 'Comparison of Rejected Transcripts {self.prefix}' --output_prefix '{self.prefix}_Rejected_Transcripts_summary_plot' --output {self.output}"
            self.process_cmd(cmd)

    def compare_parsed_summary_csv(self):
        """
        Running compare_parsed_summary_csv.py (single/paired mode)
        """
        if self.args.single:
            logging.info("Running compare_parsed_summary_csv.py (single)")
            if self.did_lifton:
                cmd = f"compare_parsed_summary_csv --single --gff_file {self.lifton_sorted_gff} --csv_file {os.path.join(self.analysis_dir, f"mikado_prepare_lifton", f"mikado_prepare_lifton_parsed_summary.csv")} --rm_prefix '1_' --output {self.output} --exclude '{self.args.exclude}'"
                self.process_cmd(cmd)
            else:
                cmd = f"compare_parsed_summary_csv --single --gff_file {self.liftoff_sorted_gff} --csv_file {os.path.join(self.analysis_dir, f"mikado_prepare_liftoff", f"mikado_prepare_liftoff_parsed_summary.csv")} --rm_prefix '1_' --output {self.output} --exclude '{self.args.exclude}'"
                self.process_cmd(cmd)
        else:
            logging.info("Running compare_parsed_summary_csv.py (paired)")
            cmd = f"compare_parsed_summary_csv --gff_files {self.lifton_sorted_gff} {self.liftoff_sorted_gff} --csv_files {os.path.join(self.analysis_dir, f"mikado_prepare_lifton", f"mikado_prepare_lifton_parsed_summary.csv")} {os.path.join(self.analysis_dir, f"mikado_prepare_liftoff", f"mikado_prepare_liftoff_parsed_summary.csv")} --rm_prefix '1_' --labels lifton,liftoff --output_prefix {self.prefix} --output {self.output} --plot_title 'Comparison of Retained and Rejected IDs {self.prefix}' --exclude '{self.args.exclude}'"
            self.process_cmd(cmd)

    def extract_gene_id_gff(self, label):
        """
        Running extract_gene_id_gff for retained IDs
        """
        id_file = os.path.join(self.analysis_dir, f"{label}.sorted.ids.retained.txt")
        out_tsv = os.path.join(self.output, f"{label}.sorted.ids.retained.tsv")
        logging.info(f"Running extract_gene_id_gff for {label}")
        cmd = None
        if label == "lifton":
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
        if label == "lifton":
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
                f" - lifton retained GFF: {os.path.join(self.output, f'lifton.sorted.retained.gff')}"
            )
            print(
                f" - lifton retained TSV: {os.path.join(self.output, f'lifton.sorted.ids.retained.tsv')}"
            )
        if self.did_liftoff:
            print(
                f" - liftoff retained GFF: {os.path.join(self.output, f'liftoff.sorted.retained.gff')}"
            )
            print(
                f" - liftoff retained TSV: {os.path.join(self.output, f'liftoff.sorted.ids.retained.tsv')}"
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
        logging.info("Starting filter-lifton-liftoff")
        logging.info(f"Command: {executed_command}")

        #######################################
        # 7) process_gff:
        #    strand-check → gffread → sort → mikado prepare
        #######################################
        if self.lifton_gff:
            logging.info("Processing Lifton GFF file")
            self.process_gff(self.lifton_gff, label="lifton")
        if self.liftoff_gff:
            logging.info("Processing Liftoff GFF file")
            self.process_gff(self.liftoff_gff, label="liftoff")

        #######################################
        # 8) Generate summary_stats.csv via parse_prepare_log
        #######################################
        self.generate_summary_stats()

        #######################################
        # 9) Prepend strand-check summary & recalc Total
        #######################################
        if self.did_lifton:
            self.append_strand_summary(label="lifton")
        if self.did_liftoff:
            self.append_strand_summary(label="liftoff")

        #######################################
        # 10) Regenerate summary plot
        #######################################
        self.regenerate_summary_plot()

        #######################################
        # 11) compare_parsed_summary_csv.py
        #######################################
        self.compare_parsed_summary_csv()

        #######################################
        # 12) extract_gene_id_gff.py
        #######################################
        if self.lifton_gff:
            self.extract_gene_id_gff(label="lifton")
        if self.liftoff_gff:
            self.extract_gene_id_gff(label="liftoff")

        #######################################
        # 13) mikado util grep
        #######################################
        if self.lifton_gff:
            self.mikado_util_grep(label="lifton")
        if self.liftoff_gff:
            self.mikado_util_grep(label="liftoff")

        #######################################
        # 14) Final summary
        #######################################
        self.final_summary()

        # Here would be the main logic of the filtering and comparison
        logging.info("Finished filter-lifton-liftoff")


def main():
    parser = argparse.ArgumentParser(
        prog=script,
        formatter_class=HelpFormatter,
        description="""
        Lifton/Liftoff transcript filtering and comparison pipeline
        """,
        epilog=f"{help_text_note}\nContact: {__author__} ({__email__})",
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
    # --prefix for plots titles
    parser.add_argument(
        "-p",
        "--prefix",
        help="Provide a label for plot titles to distinguish images from multiple runs. If not provided, defaults to base name of the --output. For example, you can use 'Wheat_Accession1' or 'SpeciesX' [default:%(default)s]",
    )
    # --output directory
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="Provide output directory [default:%(default)s]",
    )
    # --single
    parser.add_argument(
        "-s",
        "--single",
        action="store_true",
        help="Only one GFF allowed. Default is to process both --lifton_gff and --liftoff_gff inputs [default:%(default)s]",
    )
    # add threads option
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        help="Number of threads to use for copying",
    )
    # --exclude
    parser.add_argument(
        "-e",
        "--exclude",
        type=str,
        default="",
        help="Provide a comma-separated list of types to exclude from processing. Exclude categories (default: none, e.g. 'Size under minimum,Redundant' etc [default:%(default)s]",
    )
    # --minimum-cdna-length
    parser.add_argument(
        "-m",
        "--minimum-cdna-length",
        type=int,
        default=48,
        help="Provide minimum cDNA length for filtering [default:%(default)s]",
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

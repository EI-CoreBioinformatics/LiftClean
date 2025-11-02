#!/usr/bin/env python3

import re
import sys
import os
import argparse
from collections import defaultdict


def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="Check GFF files for strand inconsistencies and child-parent seqID mismatches across gene and transcript types."
    )
    parser.add_argument(
        "gff_file", metavar="gff_file", type=str, help="Path to GFF file to be checked."
    )
    parser.add_argument(
        "--remove",
        action="store_true",
        help="Remove inconsistent features and output a corrected GFF file.",
    )
    parser.add_argument(
        "-od",
        "--output_dir",
        type=str,
        default=None,
        help="Optional output directory for corrected GFF file (requires --remove).",
    )
    parser.add_argument(
        "--check_parent_seq_id",
        action="store_true",
        help="Check for child features whose sequence-ID differs from their parent.",
    )
    args = parser.parse_args()

    gff_file_path = args.gff_file
    if not os.path.isfile(gff_file_path):
        print(f"Error: File '{gff_file_path}' not found.")
        sys.exit(1)

    # Data structures
    gene_strands = defaultdict(lambda: {"+": set(), "-": set()})
    parent_strands = defaultdict(lambda: {"+": set(), "-": set()})
    gene_to_transcripts = defaultdict(list)
    transcript_children = defaultdict(list)
    feature_lines = []
    transcript_to_gene = {}
    lines_to_remove = set()
    affected_genes = set()

    # Lookup maps
    seqid_map = {}  # seqid_map[line_num] = seqID
    id_to_line = {}  # id_to_line[featureID] = line_num

    # For summary reporting
    removal_log = defaultdict(set)  # subcategory -> transcript IDs

    # Read and classify features
    with open(gff_file_path, "r") as infile:
        for line_num, raw in enumerate(infile, start=1):
            line = raw.rstrip("\n")
            feature_lines.append((line_num, line))
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 1:
                continue
            seqid_map[line_num] = parts[0]
            if len(parts) < 9:
                continue

            ftype, strand, attrs = parts[2], parts[6], parts[8]

            # Extract ID if present
            id_match = re.search(r"ID=([^;]+)", attrs)
            fid = id_match.group(1) if id_match else None
            if fid:
                id_to_line[fid] = line_num

            # Extract Parent if present
            parent_match = re.search(r"Parent=([^;]+)", attrs)
            pid = parent_match.group(1) if parent_match else None

            # Classify transcripts: any feature with an ID and a Parent, excluding exon/CDS
            if fid and pid and ftype not in ["exon", "CDS"]:
                transcript_to_gene[fid] = pid
                gene_to_transcripts[pid].append(line_num)
                if strand in ["+", "-"]:
                    gene_strands[pid][strand].add(line_num)

            # Classify exon/CDS children
            if pid and ftype in ["exon", "CDS"]:
                transcript_children[pid].append(line_num)

            # Track strand of any child under its parent
            if pid and strand in ["+", "-"]:
                parent_strands[pid][strand].add(line_num)

    # 1) Strand consistency: transcripts under genes
    for gene_id, strands in gene_strands.items():
        line_num = id_to_line.get(gene_id)
        if not line_num:
            continue
        parts = feature_lines[line_num - 1][1].split("\t")
        if len(parts) < 3 or parts[2] not in ("gene", "ncRNA_gene"):
            continue  # only check for real gene entries
        if strands["+"] and strands["-"]:
            print(f"Error: Gene '{gene_id}' has transcripts on both + and - strands.")
            for ln in sorted(strands["+"]):
                print(f"  + (line {ln}): {feature_lines[ln - 1][1]}")
            for ln in sorted(strands["-"]):
                print(f"  - (line {ln}): {feature_lines[ln - 1][1]}")
            print()
            if args.remove:
                lines_to_remove.update(strands["+"] | strands["-"])
                affected_genes.add(gene_id)
                for ln in strands["+"] | strands["-"]:
                    attrs = feature_lines[ln - 1][1].split("\t")[8]
                    tm = re.search(r"ID=([^;]+)", attrs)
                    if tm:
                        tid = tm.group(1)
                        lines_to_remove.update(transcript_children[tid])
                        removal_log["Strand conflict gene"].add(tid)

    # 2) Strand consistency: child features under transcripts (skip genes)
    for parent_id, strands in parent_strands.items():
        parent_ln = id_to_line.get(parent_id)
        if parent_ln:
            parent_type = feature_lines[parent_ln - 1][1].split("\t")[2]
            if parent_type in ("gene", "ncRNA_gene"):
                continue  # already handled above

        if strands["+"] and strands["-"]:
            print(f"Error: Parent '{parent_id}' has child features on both strands.")
            for ln in sorted(strands["+"]):
                print(f"  + (line {ln}): {feature_lines[ln - 1][1]}")
            for ln in sorted(strands["-"]):
                print(f"  - (line {ln}): {feature_lines[ln - 1][1]}")
            print()
            if args.remove:
                lines_to_remove.update(strands["+"] | strands["-"])
                lines_to_remove.add(min(strands["+"] | strands["-"]))
                if parent_id in transcript_to_gene:
                    affected_genes.add(transcript_to_gene[parent_id])
                removal_log["Strand conflict child"].add(parent_id)

    # 3) Sequence-ID mismatches
    if args.check_parent_seq_id:
        for ln, raw in feature_lines:
            if raw.startswith("#"):
                continue
            parts = raw.split("\t")
            if len(parts) < 9:
                continue
            sid, attrs = parts[0], parts[8]
            pm = re.search(r"Parent=([^;]+)", attrs)
            if not pm:
                continue
            pid = pm.group(1)
            if pid in id_to_line:
                pln = id_to_line[pid]
                ps = seqid_map.get(pln)
                if ps and ps != sid:
                    print(
                        f"Error: Child on line {ln} (seqID='{sid}') differs from parent '{pid}' on line {pln} (seqID='{ps}')."
                    )
                    print(f"  Child:  {raw}")
                    print(f"  Parent: {feature_lines[pln - 1][1]}\n")
                    if args.remove:
                        lines_to_remove.add(ln)
                        if pid in transcript_to_gene:
                            affected_genes.add(transcript_to_gene[pid])
                        removal_log["Seqid mismatch"].add(pid)

    # 4) Cascade: remove transcripts with no exon/CDS
    for tid, children in transcript_children.items():
        if all(ln in lines_to_remove for ln in children):
            print(f"All exon/CDS removed for transcript {tid}; removing transcript.")
            if tid in id_to_line:
                lines_to_remove.add(id_to_line[tid])

    # 5) Cascade: remove genes with no transcripts
    for gene_id in set(affected_genes):
        remaining = [
            ln for ln in gene_to_transcripts[gene_id] if ln not in lines_to_remove
        ]
        if not remaining:
            print(f"No remaining transcripts for gene {gene_id}; removing gene.")
            if gene_id in id_to_line:
                lines_to_remove.add(id_to_line[gene_id])

    # 6) Adjust gene coordinates based on remaining transcripts
    if args.remove:
        for gene_id in affected_genes:
            remaining = [
                ln for ln in gene_to_transcripts[gene_id] if ln not in lines_to_remove
            ]
            if remaining and gene_id in id_to_line:
                new_start = min(
                    int(feature_lines[ln - 1][1].split("\t")[3]) for ln in remaining
                )
                new_end = max(
                    int(feature_lines[ln - 1][1].split("\t")[4]) for ln in remaining
                )
                gln = id_to_line[gene_id]
                orig = feature_lines[gln - 1][1].split("\t")
                ostart, oend = int(orig[3]), int(orig[4])
                if (new_start, new_end) != (ostart, oend):
                    print(f"Updating gene {gene_id} coords to {new_start}-{new_end}")
                    orig[3], orig[4] = str(new_start), str(new_end)
                    feature_lines[gln - 1] = (gln, "\t".join(orig))
                    print(f"Updated line {gln}: {feature_lines[gln - 1][1]}")

    # Write corrected GFF
    if args.remove:
        base = os.path.basename(os.path.splitext(gff_file_path)[0]) + ".corrected.gff"
        out_path = (
            os.path.join(args.output_dir, base)
            if args.output_dir
            else f"{os.path.splitext(gff_file_path)[0]}.corrected.gff"
        )
        os.makedirs(
            os.path.dirname(out_path), exist_ok=True
        ) if args.output_dir else None
        with open(out_path, "w") as out:
            for ln, line in feature_lines:
                if ln not in lines_to_remove:
                    out.write(line + "\n")
        print(f"\nCorrected GFF file saved as '{out_path}'")

        # Summary
        if removal_log:
            print("\nSummary:")
            print("Category,SubCategory,Count,Unique TID Count")
            for subcat, tids in removal_log.items():
                print(f"discarding,{subcat},{len(tids)},{len(set(tids))}")

    print("Checks complete.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import re
import sys
import os
import argparse
from collections import defaultdict, deque


def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description=(
            "Check GFF files for: (1) strand inconsistencies; "
            "(2) child-parent seqID mismatches; and "
            "(3) duplicate IDs limited to parent features (opt-in). "
            "Optionally remove inconsistent/duplicate features and write a corrected GFF."
        )
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
    parser.add_argument(
        "--check_dup_ids",
        action="store_true",
        help="Check for duplicate IDs, but only for features that are parents; if found, remove all instances and all descendants.",
    )
    parser.add_argument(
        "--allow_trans_splicing",
        action="store_true",
        help="If set, skip all checks for any parent (and its children) where the parent has 'exception=trans-splicing'.",
    )
    args = parser.parse_args()

    gff_file_path = args.gff_file
    if not os.path.isfile(gff_file_path):
        print(f"Error: File '{gff_file_path}' not found.")
        sys.exit(1)

    # Data structures
    gene_strands = defaultdict(lambda: {"+": set(), "-": set()})
    parent_strands = defaultdict(lambda: {"+": set(), "-": set()})
    gene_to_transcripts = defaultdict(list)  # gene_id -> [transcript line nums]
    transcript_children = defaultdict(
        list
    )  # transcript_id -> [exon/CDS line nums] (legacy cascade)
    parent_to_children = defaultdict(
        list
    )  # ANY parent_id -> [child line nums] (global cascade graph)
    feature_lines = []  # [(line_num, raw_line)]
    transcript_to_gene = {}  # transcript_id -> gene_id
    lines_to_remove = set()
    affected_genes = set()

    # Lookup maps
    seqid_map = {}  # line_num -> seqID
    id_to_line = {}  # featureID -> (last) line_num
    id_occurrences = defaultdict(list)  # featureID -> [all line nums]
    line_to_id = {}  # line_num -> featureID (if present)
    line_to_type = {}  # line_num -> feature type string (gene/mRNA/exon/CDS/...)
    id_to_type = {}  # featureID -> feature type (best-effort; last one wins)
    exempt_parent_ids = (
        set()
    )  # IDs of parents to be skipped when --allow_trans_splicing is on

    # For summary reporting (transcript-style buckets)
    # NOTE: For the 'Strand conflict gene-child' case, we count GENE IDs in this transcript-style bucket per request.
    removal_log = defaultdict(set)  # subcategory -> set of "transcript-like" IDs

    # Subcategories (fixed list so we can always print zeros)
    ALL_SUBCATS = [
        "Duplicate parent feature ID",
        "Strand conflict gene",
        "Strand conflict child",
        "Seqid mismatch",
        "Strand conflict gene-child",
    ]

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
            line_to_type[line_num] = ftype

            # Extract ID if present
            id_match = re.search(r"ID=([^;]+)", attrs)
            fid = id_match.group(1) if id_match else None
            if fid:
                id_to_line[fid] = line_num
                id_occurrences[fid].append(line_num)
                line_to_id[line_num] = fid
                id_to_type[fid] = ftype

            # Track trans-splicing exemption
            if (
                args.allow_trans_splicing
                and fid
                and "exception=trans-splicing" in attrs
            ):
                # If any occurrence of an ID has exception=trans-splicing, exempt that parent ID
                exempt_parent_ids.add(fid)

            # Extract Parent if present
            parent_match = re.search(r"Parent=([^;]+)", attrs)
            pid = parent_match.group(1) if parent_match else None

            # Build global parent graph (for any parent, not just transcripts)
            if pid:
                parent_to_children[pid].append(line_num)

            # Classify transcripts: any feature with an ID and a Parent, excluding exon/CDS
            if fid and pid and ftype not in ["exon", "CDS"]:
                transcript_to_gene[fid] = pid
                gene_to_transcripts[pid].append(line_num)
                if strand in ["+", "-"]:
                    gene_strands[pid][strand].add(line_num)

            # Classify exon/CDS children per transcript
            if pid and ftype in ["exon", "CDS"]:
                transcript_children[pid].append(line_num)

            # Track strand of any child under its parent
            if pid and strand in ["+", "-"]:
                parent_strands[pid][strand].add(line_num)

    # Helpers
    def id_of_line(ln):
        return line_to_id.get(ln)

    def type_of_line(ln):
        return line_to_type.get(ln)

    def children_of_id(fid):
        return parent_to_children.get(fid, [])

    def any_exon_or_cds_child(fid):
        for cln in children_of_id(fid):
            ctype = type_of_line(cln)
            if ctype in ("exon", "CDS"):
                return True
        return False

    def is_exempt_parent(pid):
        # A parent is considered exempt if:
        #   - --allow_trans_splicing is on AND
        #   - the parent ID is in exempt_parent_ids (any of its occurrences had exception=trans-splicing)
        return args.allow_trans_splicing and pid in exempt_parent_ids

    # 0) Duplicate ID detection & cascading removals (only if enabled, and only for PARENT features)
    if args.check_dup_ids:
        # Only consider IDs that have children (i.e., appear as a Parent= somewhere)
        duplicated_parent_ids = {
            fid
            for fid, occ in id_occurrences.items()
            if len(occ) > 1 and fid in parent_to_children
        }
        # Exempt parents (and their subtrees) if requested
        duplicated_parent_ids = {
            fid for fid in duplicated_parent_ids if not is_exempt_parent(fid)
        }

        if duplicated_parent_ids:
            for dup_id in sorted(duplicated_parent_ids):
                occ = sorted(id_occurrences[dup_id])
                print(
                    f"Error: Parent feature ID '{dup_id}' appears {len(occ)} times; removing all instances and descendants."
                )
                for ln in occ:
                    print(f"  dup (line {ln}): {feature_lines[ln - 1][1]}")

                if args.remove:
                    # Remove all duplicates and their descendants (BFS/DFS over parent_to_children)
                    to_visit = deque()

                    # Seed with all occurrences
                    for ln in occ:
                        lines_to_remove.add(ln)

                    # Enqueue direct children (ID-based graph)
                    for cln in parent_to_children.get(dup_id, []):
                        to_visit.append(cln)

                    while to_visit:
                        cur = to_visit.popleft()
                        if cur in lines_to_remove:
                            continue
                        # Skip exempt children if their own ID is an exempt parent
                        cur_id = id_of_line(cur)
                        if cur_id and is_exempt_parent(cur_id):
                            continue
                        lines_to_remove.add(cur)
                        if cur_id and cur_id in parent_to_children:
                            # enqueue descendants
                            for gchild in parent_to_children[cur_id]:
                                to_visit.append(gchild)

                    # Count transcript-level removals caused by this duplicate parent ID.
                    # Only count transcript IDs that are parents of exon/CDS.
                    removed_ids = set()
                    for ln, _ in feature_lines:
                        if ln in lines_to_remove:
                            rid = id_of_line(ln)
                            if rid:
                                removed_ids.add(rid)
                    for tid in removed_ids:
                        if tid in transcript_to_gene and any_exon_or_cds_child(tid):
                            removal_log["Duplicate parent feature ID"].add(tid)

                print()

    # 1) Strand consistency: transcripts under genes
    for gene_id, strands in gene_strands.items():
        if is_exempt_parent(gene_id):
            continue
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
                        # Not exempting children here; if gene is not exempt, its transcripts count
                        lines_to_remove.update(transcript_children[tid])
                        removal_log["Strand conflict gene"].add(tid)

    # 1b) Strand consistency: direct children under genes (e.g., exon/CDS directly parented by gene)
    #     IMPORTANT: consider ALL occurrences of the gene ID (duplicates) when comparing strands.
    for gene_id in list(parent_to_children.keys()):
        if is_exempt_parent(gene_id):
            continue

        occ = id_occurrences.get(gene_id, [])
        if not occ:
            continue

        # Collect strands for ALL occurrences that are actual gene/ncRNA_gene rows
        gene_occ_strands = []
        for gln in occ:
            parts = feature_lines[gln - 1][1].split("\t")
            if len(parts) < 7:
                continue
            if parts[2] not in ("gene", "ncRNA_gene"):
                continue
            strand = parts[6]
            if strand in ["+", "-"]:
                gene_occ_strands.append(strand)

        if not gene_occ_strands:
            continue  # no valid strand-bearing gene rows to compare against

        # Examine direct children of the (logical) gene ID
        bad_children = []
        for cln in parent_to_children.get(gene_id, []):
            cparts = feature_lines[cln - 1][1].split("\t")
            if len(cparts) < 7:
                continue
            ctype, cstrand = cparts[2], cparts[6]
            if ctype not in ("exon", "CDS"):
                continue
            if cstrand not in ["+", "-"]:
                continue

            # If ANY gene occurrence has a different strand than this child, it’s a conflict
            if any(gs != cstrand for gs in gene_occ_strands):
                bad_children.append(cln)

        if bad_children:
            parent_strand_summary = ",".join(sorted(set(gene_occ_strands)))
            print(
                f"Error: Gene '{gene_id}' (strand(s) {parent_strand_summary}) has direct child features on the opposite strand."
            )
            for ln in sorted(bad_children):
                print(f"  child (line {ln}): {feature_lines[ln - 1][1]}")
            print()
            if args.remove:
                # Remove those opposite-strand children and any of their descendants
                stack = list(bad_children)
                while stack:
                    cur = stack.pop()
                    if cur in lines_to_remove:
                        continue
                    cur_id = id_of_line(cur)
                    # If descendant is an exempt parent, skip removing it and its subtree
                    if cur_id and is_exempt_parent(cur_id):
                        continue
                    lines_to_remove.add(cur)
                    if cur_id and cur_id in parent_to_children:
                        stack.extend(parent_to_children[cur_id])

                # Ensure step 5 evaluates this gene for removal
                affected_genes.add(gene_id)

                # Count genes as "transcripts" for this specific subcategory, per downstream expectations
                removal_log["Strand conflict gene-child"].add(gene_id)

    # 2) Strand consistency: child features under transcripts (skip genes)
    for parent_id, strands in parent_strands.items():
        if is_exempt_parent(parent_id):
            continue
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
                # Also drop the parent itself to avoid orphan children mix
                if parent_ln:
                    lines_to_remove.add(parent_ln)
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
            if is_exempt_parent(pid):
                continue  # skip checking children of exempt parents
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
        if not children:
            continue
        if all(ln in lines_to_remove for ln in children):
            print(f"All exon/CDS removed for transcript {tid}; removing transcript.")
            if tid in id_to_line:
                lines_to_remove.add(id_to_line[tid])

    # 4b) General orphan cleanup: if a feature's Parent was removed by any step, drop the child too.
    #     We iterate until no change to ensure deep cascades.
    if args.remove:
        changed = True
        while changed:
            changed = False
            for ln, raw in feature_lines:
                if ln in lines_to_remove:
                    continue
                if raw.startswith("#"):
                    continue
                parts = raw.split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                pm = re.search(r"Parent=([^;]+)", attrs)
                if not pm:
                    continue
                pid = pm.group(1)
                # If any occurrence of this parent ID is marked for removal (i.e., we removed the parent entity),
                # then drop the child. Exempt parents won't be marked removed, so their children stay.
                parent_removed = any(
                    oln in lines_to_remove for oln in id_occurrences.get(pid, [])
                )
                if parent_removed:
                    lines_to_remove.add(ln)
                    changed = True

    # 5) Cascade: remove genes with no transcripts
    for gene_id in set(affected_genes):
        remaining = [
            ln for ln in gene_to_transcripts[gene_id] if ln not in lines_to_remove
        ]
        if not remaining:
            print(f"No remaining transcripts for gene {gene_id}; removing gene.")
            # Remove ALL occurrences of this gene ID (not just last)
            for gln in id_occurrences.get(gene_id, []):
                lines_to_remove.add(gln)

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
        # Only make dir if user provided output_dir
        if args.output_dir:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as out:
            for ln, line in feature_lines:
                if ln not in lines_to_remove:
                    out.write(line + "\n")
        print(f"\nCorrected GFF file saved as '{out_path}'")

        # Summary (always print all subcategories with zeroes when absent)
        print("\nSummary:")
        print("Category,SubCategory,Count,Unique TID Count")
        for subcat in ALL_SUBCATS:
            ids = removal_log.get(subcat, set())
            print(f"discarding,{subcat},{len(ids)},{len(set(ids))}")

    print("Checks complete.")


if __name__ == "__main__":
    main()

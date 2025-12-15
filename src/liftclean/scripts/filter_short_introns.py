#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter GFF3 on intron size with robust gene/transcript handling.

- Gene-level types configurable via --gene_types (default: gene,ncRNA_gene,pseudogene).
- Transcript-level types configurable via --transcript_types + regex
  (default matches *RNA and 'transcript'; includes primary_transcript and pseudogenic_transcript).
- Creates a synthetic transcript for genes that have children directly under the gene:
    * If the parent gene is a pseudogene, the synthetic feature type is 'pseudogenic_transcript'.
    * Otherwise the synthetic feature type is 'transcript'.
- Re-parents single-parent children (exon/CDS/UTRs) from gene -> synthetic transcript.
- Skips transcripts/children whose parent gene type is NOT in --gene_types (prevents orphans).
- Supports transcript->transcript nesting: a transcript is kept if its *root* parent is an allowed gene.
- Logs intron statistics only over multi-exonic transcripts.
"""

import argparse
import os
import re
import logging
import sys
from collections import defaultdict
import statistics

# GFF3 column indices
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTE = range(9)

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)

# Precompiled regex for speed
ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")

# Child feature types under transcripts/genes we may need to (re)parent or skip
CHILD_FEATURE_TYPES = {
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "five_prime_utr",
    "three_prime_utr",
}


def fast_get(attr: str, rx: re.Pattern):
    """Fast attribute extractor: returns first group if match else None."""
    m = rx.search(attr)
    return m.group(1) if m else None


def get_id(file, line, attribute, field):
    """
    Extract value for 'field=' from the attributes column (verbose errors).
    Used in less-hot paths; hot-paths use fast_get().
    """
    pattern = field + "=([^;]+)"
    id_field = None
    try:
        m = re.search(pattern, attribute)
        id_field = m.group(1) if m else None
        if id_field is None:
            logging.error(
                f"Error: Cannot extract attribute type '{field}=' from the file '{file}'. "
                f"Please check if the attribute is present in line:\n{line}\n"
            )
    except AttributeError as err:
        logging.error(
            f"Error: {err}. '{field}' cannot be extracted from the file '{file}' line below\n{line}"
        )

    if id_field is None:
        logging.error(f"Error: Cannot find {field} from the file '{file}', exiting..")
        sys.exit(1)
    return id_field


def is_transcript(type_str, transcript_types, transcript_regex):
    """Return True if a feature type should be treated as a transcript."""
    if type_str in transcript_types:
        return True
    if transcript_regex and re.search(transcript_regex, type_str):
        return True
    return False


def is_gene_type(type_str, gene_types):
    """Return True if a feature type should be treated as a gene."""
    return type_str in gene_types


def get_gff3_info(filename, transcript_types, transcript_regex, gene_types):
    logging.info("Processing file {filename} ...".format(filename=filename))

    gene_info = defaultdict(
        dict
    )  # allowed genes only (per --gene_types), stores start/end/type
    mrna_info = defaultdict(dict)  # transcript-like features (per is_transcript)
    exon_info = defaultdict(list)  # parent_id -> exon coords (parent may be tx or gene)
    gene_to_transcripts = defaultdict(list)  # gene_id -> list(tids)
    children_under_gene_coords = defaultdict(list)  # gene_id -> list[[start,end],...]

    # Book-keeping for logging and validation
    id_to_type = {}  # any feature with an ID -> its type
    nonstandard_gene_type_counts = defaultdict(
        int
    )  # types that include 'gene' but not in --gene_types

    with open(filename, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            # STRICT TAB-DELIMITED parsing
            x = line.split("\t")
            if len(x) != 9:
                logging.error(
                    "Error: Not a standard GFF3 9-column tab-delimited line:\n{line}\nExiting...".format(
                        line=line
                    )
                )
                sys.exit(1)

            # Record ID->type if present (without emitting errors for features lacking ID)
            _id = fast_get(x[ATTRIBUTE], ID_RE)
            if _id:
                id_to_type[_id] = x[TYPE]

            # normalized start/end
            s_raw, e_raw = x[START], x[END]
            start, end = (
                (int(s_raw), int(e_raw))
                if int(s_raw) < int(e_raw)
                else (int(e_raw), int(s_raw))
            )

            # Gene-level features (allowed)
            if is_gene_type(x[TYPE], gene_types):
                gene_id = _id if _id else get_id(filename, line, x[ATTRIBUTE], "ID")
                if gene_id in gene_info:
                    logging.error(
                        "Error: Duplicate gene encountered '{gene_id}'. Please sort GFF3 (e.g., with genometools). Exiting...".format(
                            gene_id=gene_id
                        )
                    )
                    sys.exit(1)
                gene_info[gene_id]["start"] = start
                gene_info[gene_id]["end"] = end
                gene_info[gene_id]["type"] = x[TYPE]

            # Track non-standard gene-like types (contain 'gene' but not allowed)
            elif "gene" in x[TYPE].lower():
                nonstandard_gene_type_counts[x[TYPE]] += 1

            # Transcript-like features
            if is_transcript(x[TYPE], transcript_types, transcript_regex):
                tid = _id if _id else get_id(filename, line, x[ATTRIBUTE], "ID")
                parent_gene_or_tx = fast_get(x[ATTRIBUTE], PARENT_RE) or get_id(
                    filename, line, x[ATTRIBUTE], "Parent"
                )
                if tid in mrna_info:
                    logging.error(
                        "Error: Duplicate transcript encountered '{mrna_id}'. Please sort GFF3 (e.g., with genometools). Exiting...".format(
                            mrna_id=tid
                        )
                    )
                    sys.exit(1)
                mrna_info[tid]["parent_id"] = parent_gene_or_tx
                mrna_info[tid]["start"] = start
                mrna_info[tid]["end"] = end
                # Only map to gene if the parent is a gene (may be transcript in nested cases)
                if parent_gene_or_tx in gene_types or parent_gene_or_tx in gene_info:
                    gene_to_transcripts[parent_gene_or_tx].append(tid)

            # Exons used for intron stats (note: parent may be transcript OR gene)
            if x[TYPE] == "exon":
                parent_id = fast_get(x[ATTRIBUTE], PARENT_RE) or get_id(
                    filename, line, x[ATTRIBUTE], "Parent"
                )
                exon_info[parent_id].append([start, end])
                # track exons under a gene (we'll filter to allowed genes later)
                children_under_gene_coords[parent_id].append([start, end])

            # Other child features that might be under a gene
            if x[TYPE] in CHILD_FEATURE_TYPES and x[TYPE] != "exon":
                parent_id = fast_get(x[ATTRIBUTE], PARENT_RE) or get_id(
                    filename, line, x[ATTRIBUTE], "Parent"
                )
                children_under_gene_coords[parent_id].append([start, end])

    return (
        gene_info,
        mrna_info,
        exon_info,
        gene_to_transcripts,
        children_under_gene_coords,
        id_to_type,
        nonstandard_gene_type_counts,
    )


def synthesize_gene_children_transcripts(
    gene_info,
    mrna_info,
    exon_info,
    children_under_gene_coords,
    gene_to_transcripts,
):
    """
    Create synthetic transcripts for allowed genes that have children directly under the gene.
    Move exons under gene -> exons under synthetic transcript so introns are computed.
    Also register the synthetic in mrna_info and gene_to_transcripts.

    NOTE: We only choose the synthetic *type* when writing output (knowing the parent gene type).
    Here we just allocate an ID and register coordinates so introns can be computed if needed.
    """
    synthetic_tx_id_by_gene = {}
    synthetic_tx_span_by_gene = {}

    # reserve existing IDs to avoid collisions
    existing_ids = set(mrna_info.keys())

    for gene_id, coords in children_under_gene_coords.items():
        if gene_id not in gene_info:
            continue  # only synthesize for allowed genes
        if not coords:
            continue

        # span for synthetic transcript from all direct children under gene (exons/UTRs/CDS)
        min_s = min(c[0] for c in coords)
        max_e = max(c[1] for c in coords)

        # unique synthetic transcript ID (stable prefix 'transcript-'; type is decided at emit-time)
        base = f"transcript-{gene_id}"
        synth_id = base
        suffix = 1
        while synth_id in existing_ids:
            suffix += 1
            synth_id = f"{base}.{suffix}"
        existing_ids.add(synth_id)

        synthetic_tx_id_by_gene[gene_id] = synth_id
        synthetic_tx_span_by_gene[gene_id] = (min_s, max_e)

        # register synthetic transcript as a normal transcript (so it participates in intron calc)
        mrna_info[synth_id] = {
            "parent_id": gene_id,
            "start": min_s,
            "end": max_e,
        }
        gene_to_transcripts[gene_id].append(synth_id)

        # re-key exons under the gene to the synthetic transcript so introns get computed
        if gene_id in exon_info:
            exon_info[synth_id].extend(exon_info[gene_id])
            # keep original key to be harmless; intron computation will only consider mrna_info keys

    return synthetic_tx_id_by_gene, synthetic_tx_span_by_gene


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
        if len(coords) == 1:  # mono-exonic models
            intron_info[transcript] = []
        else:
            for i, _ in enumerate(coords[:-1]):
                intron_info[transcript].append([coords[i][1] + 1, coords[i + 1][0] - 1])

    # Build stats ONLY over multi-exonic transcripts (eligible for filtering)
    intron_stats_info = defaultdict(dict)
    multi_max_list = []
    starts_multi = []
    ends_multi = []

    for transcript, coords in intron_info.items():
        coords.sort(key=lambda interval: interval[0])
        if len(coords) == 0:
            continue  # skip mono-exonic for stats
        sizes = [c[1] - c[0] + 1 for c in coords]
        intron_stats_info[transcript]["start"] = coords[0][1] - coords[0][0] + 1
        intron_stats_info[transcript]["end"] = coords[-1][1] - coords[-1][0] + 1
        intron_stats_info[transcript]["min"] = min(sizes)
        intron_stats_info[transcript]["max"] = max(sizes)

        multi_max_list.append(intron_stats_info[transcript]["max"])
        starts_multi.append(intron_stats_info[transcript]["start"])
        ends_multi.append(intron_stats_info[transcript]["end"])

    if multi_max_list:
        logging.info(
            f"# - Minimum intron size:\t{min(intron_stats_info[tid]['min'] for tid in intron_stats_info)}"
        )
        logging.info(f"# - Maximum intron size:\t{max(multi_max_list)}")
        logging.info(f"# - Mean intron size:\t{statistics.mean(multi_max_list)}")
        logging.info(f"# - Median intron size:\t{statistics.median(multi_max_list)}")
        both_ends = starts_multi + ends_multi
        logging.info(f"# - Minimum terminal intron size:\t{min(both_ends)}")
        logging.info(f"# - Maximum terminal intron size:\t{max(both_ends)}")
    else:
        logging.info("# - No multi-exonic transcripts found for intron statistics.")

    return intron_info, intron_stats_info


def root_gene_of_transcript(tid, mrna_map, gene_map, _cache=None):
    """
    Follow Parent links up from tid until we reach a gene in gene_map.
    Returns gene_id or None. Caches results for speed.
    """
    if _cache is None:
        _cache = {}
    if tid in _cache:
        return _cache[tid]
    seen = set()
    cur = tid
    while True:
        if cur in _cache:
            _cache[tid] = _cache[cur]
            return _cache[cur]
        if cur in seen:
            _cache[tid] = None  # cycle guard
            return None
        seen.add(cur)
        parent = mrna_map.get(cur, {}).get("parent_id")
        if parent is None:
            _cache[tid] = None
            return None
        if parent in gene_map:
            _cache[tid] = parent
            return parent
        # parent is another transcript — walk up
        cur = parent


def filter_intron(
    gene_info,
    mrna_info,
    exon_info,
    intron_stats_info,
    intron_size,
    intron_size_ends,
    detailed_discard_list_file,
):
    # Informative count: eligible multi-exonic transcripts
    eligible = sum(
        1 for tid, exs in exon_info.items() if len(exs) >= 2 and tid in mrna_info
    )
    logging.info(
        f"# Transcripts eligible for intron filtering (multi-exonic):\t{eligible}"
    )

    discard_list = {}
    for transcript in intron_stats_info:
        # only remove multi-exonic models (mono-exonic never added to intron_stats_info)
        if transcript not in mrna_info:
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

    logging.info(
        "Identified {count} transcripts (of total {mrna_count} transcripts / {gene_count} genes) to be removed".format(
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

    # Detailed discard list
    with open(detailed_discard_list_file, "w") as f:
        f.write(
            "\t".join(
                [
                    "transcript_id",
                    "gene_id",
                    "exon_coords",
                    "min_intron_size",
                    "max_intron_size",
                ]
            )
            + "\n"
        )
        for tid in discard_list.keys():
            f.write(
                "\t".join(
                    [
                        tid,
                        discard_list[tid],
                        str(exon_info[tid]),
                        str(intron_stats_info[tid]["min"]),
                        str(intron_stats_info[tid]["max"]),
                    ]
                )
                + "\n"
            )

    filtered_mrna_info = mrna_info.copy()
    for tid in discard_list:
        if tid in filtered_mrna_info:
            del filtered_mrna_info[tid]

    logging.info("Filter failed intron models and readjust gene coordinates ...")
    updated_gene_info = defaultdict(dict)
    for tid, info in filtered_mrna_info.items():
        start = info["start"]
        end = info["end"]
        parent_gene = info["parent_id"]
        # Only propagate coordinates to *real* genes
        if parent_gene not in gene_info:
            continue
        if parent_gene not in updated_gene_info:
            updated_gene_info[parent_gene]["start"] = start
            updated_gene_info[parent_gene]["end"] = end
        else:
            if int(start) < int(updated_gene_info[parent_gene]["start"]):
                updated_gene_info[parent_gene]["start"] = start
            if int(end) > int(updated_gene_info[parent_gene]["end"]):
                updated_gene_info[parent_gene]["end"] = end

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
    transcript_types,
    transcript_regex,
    gene_to_transcripts,
    children_under_gene_coords,
    gene_types,
    synthetic_tx_id_by_gene,
    synthetic_tx_span_by_gene,
):
    # Log how many synthetic transcripts were created
    logging.info(f"# Synthetic transcripts created:\t{len(synthetic_tx_id_by_gene)}")

    # Determine which transcripts are actually writable: a transcript is kept if its root gene is allowed
    allowed_gene_ids = set(gene_info.keys())
    root_cache = {}

    def root_gene(t):
        return root_gene_of_transcript(
            t, filtered_mrna_info, gene_info, _cache=root_cache
        )

    kept_transcript_ids = {
        t for t in filtered_mrna_info if (rg := root_gene(t)) in allowed_gene_ids
    }

    # Counters for logging
    skipped_transcripts_due_to_invalid_gene = len(filtered_mrna_info) - len(
        kept_transcript_ids
    )
    skipped_child_features_due_to_invalid_parent = 0
    reparented_children = 0
    multi_parent_children_left_unchanged = 0

    # Hoist lookups for speed
    is_kept_tx = kept_transcript_ids.__contains__
    is_allowed_gene = allowed_gene_ids.__contains__
    get_synth = synthetic_tx_id_by_gene.get
    get_synth_span = synthetic_tx_span_by_gene.get
    gene_type_set = set(gene_types)

    # write output GFF3
    with open(output_gff, "w") as f:
        f.write("##gff-version 3\n")
        skip_first = False

        with open(filename, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue

                # STRICT TAB-DELIMITED parsing in writer too
                x = line.split("\t")
                if len(x) != 9:
                    continue  # skip malformed lines

                if source:
                    x[SOURCE] = source

                feature_type = x[TYPE]

                # Genes (allowed types only)
                if feature_type in gene_type_set:
                    gene_id = fast_get(x[ATTRIBUTE], ID_RE) or get_id(
                        filename, line, x[ATTRIBUTE], "ID"
                    )
                    if is_allowed_gene(gene_id) and gene_id in updated_gene_info:
                        if skip_first:
                            f.write("###\n")
                        skip_first = True
                        new_start = updated_gene_info[gene_id]["start"]
                        new_end = updated_gene_info[gene_id]["end"]
                        fline = [
                            x[SEQID],
                            x[SOURCE],
                            feature_type,
                            str(new_start),
                            str(new_end),
                            x[SCORE],
                            x[STRAND],
                            x[PHASE],
                            x[ATTRIBUTE],
                        ]
                        f.write("\t".join(fline) + "\n")

                        # Emit precomputed synthetic transcript after gene, if any (and kept)
                        synth_id = get_synth(gene_id)
                        if synth_id and synth_id in kept_transcript_ids:
                            tstart, tend = get_synth_span(gene_id)
                            # choose synthetic feature type based on the gene type
                            gene_type = gene_info.get(gene_id, {}).get("type", "gene")
                            synth_feature_type = (
                                "pseudogenic_transcript"
                                if gene_type == "pseudogene"
                                else "transcript"
                            )
                            attrs = f"ID={synth_id};Parent={gene_id}"
                            tline = [
                                x[SEQID],
                                x[SOURCE],
                                synth_feature_type,
                                str(tstart),
                                str(tend),
                                ".",
                                x[STRAND],
                                ".",
                                attrs,
                            ]
                            f.write("\t".join(tline) + "\n")
                    continue  # handled gene line

                # Transcripts
                if is_transcript(feature_type, transcript_types, transcript_regex):
                    tid = fast_get(x[ATTRIBUTE], ID_RE) or get_id(
                        filename, line, x[ATTRIBUTE], "ID"
                    )
                    if is_kept_tx(tid):
                        f.write("\t".join(x) + "\n")
                        if add_intron:
                            introns = intron_info.get(tid, [])
                            if introns:
                                for i, coord in enumerate(introns):
                                    intronLine = [
                                        x[SEQID],
                                        x[SOURCE],
                                        "intron",
                                        str(coord[0]),
                                        str(coord[1]),
                                        ".",
                                        x[STRAND],
                                        ".",
                                        f"ID={tid}.intron{i + 1};Parent={tid}",
                                    ]
                                    f.write("\t".join(intronLine) + "\n")
                    continue

                # Child features (exon/CDS/UTRs)
                if feature_type in CHILD_FEATURE_TYPES:
                    parent_field = fast_get(x[ATTRIBUTE], PARENT_RE) or get_id(
                        filename, line, x[ATTRIBUTE], "Parent"
                    )

                    # Multi-parent? leave unchanged (as per preference) and count
                    if parent_field and "," in parent_field:
                        multi_parent_children_left_unchanged += 1
                        f.write("\t".join(x) + "\n")
                        continue

                    parent_id = parent_field
                    if parent_id and is_kept_tx(parent_id):
                        # child under a kept transcript
                        f.write("\t".join(x) + "\n")
                    elif parent_id and is_allowed_gene(parent_id):
                        synth = get_synth(parent_id)
                        if synth and synth in kept_transcript_ids:
                            # re-parent child under synthetic transcript for this gene (single replace)
                            x_copy = x[:]  # shallow copy
                            x_attr = x_copy[ATTRIBUTE]
                            x_copy[ATTRIBUTE] = x_attr.replace(
                                f"Parent={parent_id}", f"Parent={synth}", 1
                            )
                            f.write("\t".join(x_copy) + "\n")
                            reparented_children += 1
                        else:
                            # allowed gene but no (kept) synthetic — keep as-is
                            f.write("\t".join(x) + "\n")

                    else:
                        # Parent isn't a kept transcript nor an allowed gene -> skip to avoid orphaning
                        skipped_child_features_due_to_invalid_parent += 1
                    continue

                # Any other feature types: pass through unchanged
                f.write("\t".join(x) + "\n")

        f.write("###\n")

    # write discard list (headers say transcript_id)
    with open(discard_list_mapping_file, "w") as f:
        f.write("transcript_id\tgene_id\n")
        for tid in discard_list.keys():
            f.write(f"{tid}\t{discard_list[tid]}\n")

    # Append summary lines to the detailed discard list file (for compatibility with original script)
    with open(detailed_discard_list_file, "a") as f:
        f.write("\n\n")
        f.write("Category,SubCategory,Count,Unique TID Count\n")
        f.write(f"discarding,Short intron,{len(discard_list)},{len(discard_list)}\n")
        f.write("Checks complete.\n")

    # Summary logging
    logging.info("SUMMARY")
    logging.info(f"# Input file:\t{filename}")
    logging.info(f"# Input genes:\t{len(gene_info.keys())}")
    logging.info(f"# Input transcripts (including synthetic):\t{len(mrna_info.keys())}")
    logging.info(
        f"# transcripts removed based on intron filter:\t{len(discard_list.keys())}"
    )
    logging.info(
        f"# genes removed based on intron filter:\t{len(set(discard_list.values()))}"
    )

    output_gene_ids = [g for g in updated_gene_info.keys() if g in gene_info]
    logging.info(
        f"# Output genes:\t{len(output_gene_ids)} "
        f"[{(100 * len(output_gene_ids) / max(1, len(gene_info.keys()))):.2f}%]"
    )

    logging.info(
        f"# Output transcripts:\t{len(kept_transcript_ids)} "
        f"[{(100 * len(kept_transcript_ids) / max(1, len(mrna_info.keys()))):.2f}%]"
    )
    logging.info(
        f"# Child features re-parented to synthetic transcripts:\t{reparented_children}"
    )
    if multi_parent_children_left_unchanged > 0:
        logging.info(
            f"# Multi-parent child features left unchanged:\t{multi_parent_children_left_unchanged}"
        )

    skipped_transcripts_due_to_invalid_gene = len(filtered_mrna_info) - len(
        kept_transcript_ids
    )
    if skipped_transcripts_due_to_invalid_gene > 0:
        logging.warning(
            f"# Transcripts skipped due to non-standard or missing parent gene type: {skipped_transcripts_due_to_invalid_gene}"
        )
    if skipped_child_features_due_to_invalid_parent > 0:
        logging.warning(
            f"# Child features skipped due to invalid parent (non-standard gene or filtered transcript): {skipped_child_features_due_to_invalid_parent}"
        )

    logging.info(
        f"# Discarded models mapping file written to file: {discard_list_mapping_file}"
    )
    logging.info(
        f"# Discarded models detailed information written to file: {detailed_discard_list_file}"
    )
    logging.info(f"# Filtered output GFF3 written to file: {output_gff}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter GFF3 on intron size. Structure: [gene-like]->[transcript-like]->[exon/CDS/UTR]. "
            "Gene-level types are configurable (default: gene,ncRNA_gene,pseudogene). "
            "Also synthesizes a transcript feature for children attached directly to gene "
            "(uses 'pseudogenic_transcript' when parent is pseudogene)."
        )
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
        "--add_intron", action="store_true", help="Add intron type to the output GFF3"
    )
    # discard list file
    parser.add_argument(
        "--discard_list_mapping_file",
        required=True,
        help="Output TSV mapping transcript_id to gene_id for qualifying introns [required]",
    )
    # detailed discard list file
    parser.add_argument(
        "--detailed_discard_list_file",
        required=True,
        help="Output TSV listing introns meeting the filter criteria [required]",
    )
    # output file
    parser.add_argument(
        "-o", "--output", required=True, help="Output GFF3 file [required]"
    )

    # Configurable transcript and gene types
    parser.add_argument(
        "--transcript_types",
        default="mRNA,primary_transcript,transcript,lnc_RNA,ncRNA,miRNA,rRNA,tRNA,snoRNA,snRNA,scaRNA,pseudogenic_transcript,antisense_RNA",
        help="Comma-separated list of feature types to treat as transcripts (default: %(default)s)",
    )
    parser.add_argument(
        "--transcript_types_regex",
        default="(?:^|_)?RNA$|^transcript$",
        help="Regex to additionally match transcript-like types. Use '' to disable. (default: %(default)s)",
    )
    parser.add_argument(
        "--gene_types",
        default="gene,ncRNA_gene,pseudogene",
        help="Comma-separated list of feature types to treat as gene-level (default: %(default)s)",
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
    transcript_types = set([t for t in args.transcript_types.split(",") if t])
    transcript_regex = (
        args.transcript_types_regex if args.transcript_types_regex else None
    )
    gene_types = set([t for t in args.gene_types.split(",") if t])

    # create parent dir if not exists
    os.makedirs(os.path.dirname(discard_list_mapping_file), exist_ok=True)
    os.makedirs(os.path.dirname(detailed_discard_list_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_gff), exist_ok=True)

    (
        gene_info,
        mrna_info,
        exon_info,
        gene_to_transcripts,
        children_under_gene_coords,
        id_to_type,
        nonstandard_gene_type_counts,
    ) = get_gff3_info(filename, transcript_types, transcript_regex, gene_types)

    # Informative logging about non-standard gene-like types
    if nonstandard_gene_type_counts:
        total_nonstd = sum(nonstandard_gene_type_counts.values())
        details = ", ".join(
            [f"{t}:{c}" for t, c in sorted(nonstandard_gene_type_counts.items())]
        )
        logging.warning(
            f"# Found {total_nonstd} gene-like features whose type is not in --gene_types; "
            f"their descendants will be dropped. Types: {details}"
        )

    # Create synthetic transcripts BEFORE intron computation so they are filterable
    synthetic_tx_id_by_gene, synthetic_tx_span_by_gene = (
        synthesize_gene_children_transcripts(
            gene_info,
            mrna_info,
            exon_info,
            children_under_gene_coords,
            gene_to_transcripts,
        )
    )

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
        transcript_types,
        transcript_regex,
        gene_to_transcripts,
        children_under_gene_coords,
        gene_types,
        synthetic_tx_id_by_gene,
        synthetic_tx_span_by_gene,
    )


if __name__ == "__main__":
    main()

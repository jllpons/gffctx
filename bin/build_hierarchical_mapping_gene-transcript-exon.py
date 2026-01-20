#!/usr/bin/env python3
"""
RULES:
- Every output row MUST contain a gene_id.
- If a transcript has no valid gene_id mapping => omit it entirely.
- If an exon has no valid transcript mapping => omit it entirely.
- Output ONLY complete triplets: gene_id, transcript_id, exon_id
  (No blank transcript/exon columns, no "orphan" output, no inference from exon IDs.)
"""

import argparse
import csv
import sys
from collections import defaultdict


def is_blank_row(row) -> bool:
    """True if the row is empty or only contains whitespace/empty cells."""
    return (not row) or all((c is None or str(c).strip() == "") for c in row)


def clean(cell) -> str:
    """Convert a CSV cell to a stripped string (safe for None)."""
    return (cell or "").strip()


def load_gene_ids(path: str) -> set[str]:
    """
    Read --geneids CSV and return a set of valid gene IDs.
    Assumes gene ID is in column 0.
    """
    gene_ids: set[str] = set()

    with open(path, newline="") as f:
        r = csv.reader(f, delimiter=",")
        for row in r:
            if is_blank_row(row):
                continue
            gene_id = clean(row[0])
            if gene_id:
                gene_ids.add(gene_id)

    return gene_ids


def load_gene_transcript_map(path: str, valid_genes: set[str]) -> dict[str, str]:
    """
    Read --geneid_transcriptids CSV and build transcript -> gene mapping.

    Rules:
    - If gene_id is blank => ignore (transcript is "orphan", but we do NOT output orphans).
    - If transcript_id is blank => ignore.
    - If valid_genes is provided, we only accept mappings where gene_id is in valid_genes.
      (This prevents outputting unknown genes.)
    - If a transcript maps to multiple different genes, we drop it (ambiguous).
    """
    # Collect ALL gene assignments per transcript first (to detect conflicts).
    transcript_to_genes: dict[str, set[str]] = defaultdict(set)

    with open(path, newline="") as f:
        r = csv.reader(f, delimiter=",")
        for row in r:
            if is_blank_row(row):
                continue

            # Expect: gene_id, transcript_id
            if len(row) < 2:
                continue

            gene_id = clean(row[0])
            transcript_id = clean(row[1])

            if not gene_id or not transcript_id:
                continue

            # Enforce that gene_id must be a known gene from --geneids
            if valid_genes and gene_id not in valid_genes:
                continue

            transcript_to_genes[transcript_id].add(gene_id)

    # Now resolve to a strict transcript -> gene mapping, dropping ambiguous transcripts.
    transcript_to_gene: dict[str, str] = {}
    for transcript_id, genes in transcript_to_genes.items():
        if len(genes) == 1:
            transcript_to_gene[transcript_id] = next(iter(genes))
        else:
            # Ambiguous transcript mapping: transcript appears with multiple genes.
            # We omit it entirely to preserve correctness.
            print(
                f"WARNING: transcript '{transcript_id}' maps to multiple genes {sorted(genes)}; omitting transcript",
                file=sys.stderr,
            )

    return transcript_to_gene


def load_transcript_exon_map(path: str, transcript_to_gene: dict[str, str]) -> dict[str, set[str]]:
    """
    Read --transcriptid_exonids CSV and build transcript -> set(exon_id).

    Rules:
    - Only accept exon mappings where transcript_id is a transcript we already know maps to a gene.
      (If transcript isn't in transcript_to_gene, exon is ignored.)
    - If exon_id is blank => ignore.
    - If transcript_id is blank => ignore.
    """
    transcript_to_exons: dict[str, set[str]] = defaultdict(set)

    with open(path, newline="") as f:
        r = csv.reader(f, delimiter=",")
        for row in r:
            if is_blank_row(row):
                continue

            # Expect: transcript_id, exon_id
            if len(row) < 2:
                continue

            transcript_id = clean(row[0])
            exon_id = clean(row[1])

            if not transcript_id or not exon_id:
                continue

            # This is the key rule: exon must map to a transcript that maps to a gene.
            if transcript_id not in transcript_to_gene:
                continue

            transcript_to_exons[transcript_id].add(exon_id)

    return transcript_to_exons


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--geneids", required=True, help="CSV with gene IDs in column 0")
    p.add_argument("--geneid_transcriptids", required=True, help="CSV: gene_id,transcript_id")
    p.add_argument("--transcriptid_exonids", required=True, help="CSV: transcript_id,exon_id")
    args = p.parse_args()

    # 1) Load known genes (used to validate transcript->gene mappings)
    gene_ids = load_gene_ids(args.geneids)

    # 2) Build transcript -> gene mapping (strict; orphans ignored; ambiguous dropped)
    transcript_to_gene = load_gene_transcript_map(args.geneid_transcriptids, gene_ids)

    # 3) Build transcript -> exons mapping, but ONLY for transcripts that map to genes
    transcript_to_exons = load_transcript_exon_map(args.transcriptid_exonids, transcript_to_gene)

    # 4) Output ONLY complete (gene, transcript, exon) rows
    out = csv.writer(sys.stdout, delimiter=",", lineterminator="\n")

    rows = []
    for transcript_id, exons in transcript_to_exons.items():
        gene_id = transcript_to_gene[transcript_id]  # always exists by construction
        for exon_id in exons:
            rows.append((gene_id, transcript_id, exon_id))

    for gene_id, transcript_id, exon_id in sorted(rows):
        out.writerow([gene_id, transcript_id, exon_id])


if __name__ == "__main__":
    main()

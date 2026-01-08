#!/usr/bin/env python3

import argparse
from collections import defaultdict
import csv
import sys

def is_blank_row(row):
    """
    Check if a row is blank (i.e., all cells are None or empty strings).

    Examples:
    >>> is_blank_row([None, "", "   "])
    True
    >>> is_blank_row([None, "data", "   "])
    False
    """
    return all((c is None or str(c).strip() == "") for c in row)

def main():

    p = argparse.ArgumentParser(description="Left-outer-join three CSV files to build a hierarchical mapping of gene -> transcripts -> exons.")
    p.add_argument("--geneids", type=str, required=True)
    p.add_argument("--geneid_transcriptids", type=str, required=True)
    p.add_argument("--transcriptid_exonids", type=str, required=True)
    args = p.parse_args()

    geneid_path, geneid_transcriptid_path, transcriptid_exonid_path = args.geneids, args.geneid_transcriptids, args.transcriptid_exonids

    # Obtain the gene IDs from the provided file
    gene_ids = set()
    with open(geneid_path) as geneid_file:
        reader = csv.reader(geneid_file, delimiter=',')
        for row in reader:
            if not row or is_blank_row(row):
                continue
            gene_id = row[0]
            gene_ids.add(gene_id)

    # Build `gene -> [transcript1, transcript2, ...]` mapping
    gene_to_transcripts = defaultdict(list)
    orphan_transcripts = set()
    with open(geneid_transcriptid_path) as geneid_transcriptid_file:
        reader = csv.reader(geneid_transcriptid_file, delimiter=',')
        for row in reader:
            if not row or is_blank_row(row):
                continue
            gene_id = row[0]
            transcript_id = row[1]
            if not gene_id:
                orphan_transcripts.add(transcript_id)
                continue
            gene_to_transcripts[gene_id].append(transcript_id)

    # Build `transcript -> [exon1, exon2, ...]` mapping
    transcript_to_exons = defaultdict(list)
    orphan_exons = set()
    with open(transcriptid_exonid_path) as transcriptid_exonid_file:
        reader = csv.reader(transcriptid_exonid_file, delimiter=',')
        for row in reader:
            if not row or is_blank_row(row):
                continue
            transcript_id = row[0]
            exon_id = row[1]
            if not transcript_id:
                orphan_exons.add(exon_id)
                continue
            transcript_to_exons[transcript_id].append(exon_id)


    # Build the hierarchical mapping: gene -> transcripts -> exons
    hierarchical_mapping = {}
    for gene_id in gene_ids:
        transcripts = gene_to_transcripts.get(gene_id, [])
        hierarchical_mapping[gene_id] = {}
        for transcript_id in transcripts:
            exons = transcript_to_exons.get(transcript_id, [])
            hierarchical_mapping[gene_id][transcript_id] = exons

    # Output the hierarchical mapping as CSV
    out = csv.writer(sys.stdout, delimiter=',', lineterminator='\n')

    for gene_id, transcripts in hierarchical_mapping.items():
        for transcript_id, exons in transcripts.items():
            if exons:
                for exon_id in exons:
                    out.writerow([gene_id, transcript_id, exon_id])
            else:
                out.writerow([gene_id, transcript_id, ""])
        if not transcripts:
            out.writerow([gene_id, "", ""])

    for orphan_transcript_id in orphan_transcripts:
        exons = transcript_to_exons.get(orphan_transcript_id, [])
        if exons:
            for exon_id in exons:
                out.writerow(["", orphan_transcript_id, exon_id])
        else:
            out.writerow(["", orphan_transcript_id, ""])


if __name__ == "__main__":
    main()

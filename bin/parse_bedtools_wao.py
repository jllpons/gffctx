#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from urllib.parse import unquote
from typing import Optional, Generator
from enum import Enum
import csv
import sys
import os


GENE_CENTRIC_CSV = "gene_summary.csv"
EXONIC_INTERVALS_GFF3 = "exonic_intervals.gff3"
INTRONIC_INTERVALS_GFF3 = "intronic_intervals.gff3"
ANTISENSE_INTERVALS_GFF3 = "antisense_intervals.gff3"
INTERGENIC_INTERVALS_GFF3 = "intergenic_intervals.gff3"
AMBIGUOUS_INTERVALS_GFF3 = "ambiguous_intervals.gff3"

@dataclass
class Interval:
    seqid: str
    source: str
    _type: str
    start: int          # GFF 1-based inclusive — keep as-is, don't convert
    end: int
    score: str          # '.' or a number as a string
    strand: str         # '+', '-', or '.'
    phase: str         # '0', '1', '2', or '.'
    attrs: dict[str, str]

@dataclass
class WaoRow:
    hit: Interval
    feature: Interval | None
    overlap_bp: int

@dataclass
class GeneOverlap:
    gene_id: str
    gene_strand: str
    overlap_bp: int

class HitClass(Enum):
    AMBIGUOUS = "ambiguous"
    ANTISENSE = "antisense"
    EXONIC = "exonic"
    INTERGENIC = "intergenic"
    INTRONIC = "intronic"

@dataclass
class HitRecord:
    hit: Interval
    gene_overlaps: list[GeneOverlap]
    n_exon_overlap: int = 0
    hit_id: str | None = None
    hit_class: HitClass | None = None

@dataclass
class GeneSummary:
    gene_id: str
    exonic_hits: list[str]
    intronic_hits: list[str]
    antisense_hits: list[str]
    ambiguous_hits: list[str]


CLASS_TO_GFF3 = {
    HitClass.EXONIC:     EXONIC_INTERVALS_GFF3,
    HitClass.INTRONIC:   INTRONIC_INTERVALS_GFF3,
    HitClass.ANTISENSE:  ANTISENSE_INTERVALS_GFF3,
    HitClass.INTERGENIC: INTERGENIC_INTERVALS_GFF3,
    HitClass.AMBIGUOUS:  AMBIGUOUS_INTERVALS_GFF3,
}


def mk_hit_id(hit: Interval) -> str:
    """
    Make a string ID for a hit, e.g. "chr1:100-200(+)". This is not guaranteed to be unique, but should be good enough for this use case.
    """
    return f"{hit.seqid}:{hit.start}-{hit.end}({hit.strand})"


def parse_attrs(attr_string: str) -> dict[str, str]:
    """
    Parse 'ID=foo;Name=bar' into a dict.
    Handles trailing ';', URL-escaped values (per GFF3 spec),
    and lines where attr_string == '.' (returns empty dict).
    """
    if attr_string == "." or not attr_string:
        return {}
    out = {}
    for pair in attr_string.rstrip(";").split(";"):
        if "=" not in pair:
            continue
        k, v = pair.split("=", 1)
        out[k.strip()] = unquote(v.strip())
    return out


def parse_wao_row(line: str) -> WaoRow:
    """
    Split a tab-separated WAO line into a typed record.
    The hit side is always present; the feature side is None when overlap_bp == 0.
    """
    f = line.rstrip("\n").split("\t")
    hit = Interval(
            seqid=f[0],
            source=f[1],
            _type=f[2],
            start=int(f[3]),
            end=int(f[4]),
            score=f[5],
            strand=f[6],
            phase=f[7],
            attrs=parse_attrs(f[8])
            )
    overlap_bp = int(f[-1])
    if overlap_bp == 0:
        return WaoRow(hit=hit, feature=None, overlap_bp=0)
    feat = Interval(
            seqid=f[9],
            source=f[10],
            _type=f[11],
            start=int(f[12]),
            end=int(f[13]),
            score=f[14],
            strand=f[15],
            phase=f[16],
            attrs=parse_attrs(f[17])
            )
    return WaoRow(hit=hit, feature=feat, overlap_bp=overlap_bp)


def _ensure_record(hits: dict[str, HitRecord], waorow, hit_id: str) -> HitRecord:
    """
    Ensure that hits[hit_id] exists, creating a new HitRecord if necessary. Return the HitRecord.
    Thanks to this helper, we can always assume that hits[hit_id] exists
    when processing a wao row, and just update it with gene/exon info as needed.
    """
    if hit_id not in hits:
        hits[hit_id] = HitRecord(
            hit=waorow.hit,
            gene_overlaps=[],
            n_exon_overlap=0,
            hit_id=hit_id,
        )
    return hits[hit_id]

def collect_hits(genes_wao_path, exons_wao_path) -> dict[str, HitRecord]:
    """
    Collect hit records from two Bedtools WAO outputs: one for loci (genes),
    one for segments (exons).
    The keys of the returned dict are hit IDs (e.g. "chr1:100-200(+)"),
    and the values are HitRecord objects containing the hit interval,
    gene overlaps, and exon overlap count.
    """
    hits: dict[str, HitRecord] = {}

    with open(genes_wao_path) as f:
        for line in f:
            waorow = parse_wao_row(line)
            record = _ensure_record(hits, waorow, mk_hit_id(waorow.hit))
            if waorow.feature is not None:
                record.gene_overlaps.append(GeneOverlap(
                    gene_id=waorow.feature.attrs["ID"],
                    gene_strand=waorow.feature.strand,
                    overlap_bp=waorow.overlap_bp,
                ))

    with open(exons_wao_path) as f:
        for line in f:
            waorow = parse_wao_row(line)
            record = _ensure_record(hits, waorow, mk_hit_id(waorow.hit))
            if waorow.feature is not None:
                record.n_exon_overlap += 1

    return hits

def classify_hit(record: HitRecord) -> HitClass:
    """
    Classify a hit as one of Ambiguous, Antisense, Exonic, Intronic, or Intergenic.
    """

    if not record.gene_overlaps:
        return HitClass.INTERGENIC

    if len(record.gene_overlaps) > 1:
        return HitClass.AMBIGUOUS

    gene_overlap = record.gene_overlaps[0]
    if (record.hit.strand == '.' or gene_overlap.gene_strand == '.') or (record.hit.strand == gene_overlap.gene_strand):
        if record.n_exon_overlap > 0:
            return HitClass.EXONIC
        else:
            return HitClass.INTRONIC
    else:
        return HitClass.ANTISENSE


def interval_to_gff3_str(interval: Interval) -> str:
    """
    Convert an Interval to a GFF3 line (without the final newline).
    """
    attr_str = ";".join(f"{k}={v}" for k, v in interval.attrs.items())
    return "\t".join([
        interval.seqid,
        interval.source,
        interval._type,
        str(interval.start),
        str(interval.end),
        interval.score,
        interval.strand,
        interval.phase,
        attr_str
    ])


def assign_hit_classes(hits: dict[str, HitRecord]) -> None:
    """
    Walk through the hits and assign a hit_class to each HitRecord in-place.
    This modifies the hits dict, but does not return anything.
    """
    for record in hits.values():
        record.hit_class = classify_hit(record)


def _ensure_gene_summary(gene_summaries: dict[str, GeneSummary], gene_id: str) -> GeneSummary:
    """
    Ensure that gene_summaries[gene_id] exists, creating a new GeneSummary if necessary.
    Return the GeneSummary.
    """
    if gene_id not in gene_summaries:
        gene_summaries[gene_id] = GeneSummary(
            gene_id=gene_id,
            exonic_hits=[],
            intronic_hits=[],
            antisense_hits=[],
            ambiguous_hits=[],
        )
    return gene_summaries[gene_id]



def collect_gene_summaries(hits: dict[str, HitRecord]) -> dict[str, GeneSummary]:
    """
    Walk through the hits and collect gene-centric summaries. 
    The keys of the returned dict are gene IDs, and the values are GeneSummary objects 
    containing the counts of exonic, intronic, antisense, and ambiguous hits for each gene.
    """

    gene_summaries = {}

    for hit_id, record in hits.items():

        for gene_overlap in record.gene_overlaps:
            gene_summary = _ensure_gene_summary(gene_summaries, gene_overlap.gene_id)

            match record.hit_class:

                case HitClass.EXONIC:
                    gene_summary.exonic_hits.append(hit_id)
                case HitClass.INTRONIC:
                    gene_summary.intronic_hits.append(hit_id)
                case HitClass.ANTISENSE:
                    gene_summary.antisense_hits.append(hit_id)
                case HitClass.AMBIGUOUS:
                    gene_summary.ambiguous_hits.append(hit_id)
                case HitClass.INTERGENIC:
                    pass
                case _:
                    raise ValueError(f"Unexpected hit class {record.hit_class} for hit {hit_id}")

    return gene_summaries


def gene_summary_as_csv(gene_summaries: dict[str, GeneSummary], file=f"{os.getcwd()}/{GENE_CENTRIC_CSV}") -> None:
    """
    """
    with open(file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "gene_id",
            "n_exonic_hits",
            "n_intronic_hits",
            "n_antisense_hits",
            "n_ambiguous_hits",
            "exonic_hit_ids",
            "intronic_hit_ids",
            "antisense_hit_ids",
            "ambiguous_hit_ids",
            ])
        for gene_summary in gene_summaries.values():
            writer.writerow([
                gene_summary.gene_id,
                len(gene_summary.exonic_hits),
                len(gene_summary.intronic_hits),
                len(gene_summary.antisense_hits),
                len(gene_summary.ambiguous_hits),
                ";".join(gene_summary.exonic_hits),
                ";".join(gene_summary.intronic_hits),
                ";".join(gene_summary.antisense_hits),
                ";".join(gene_summary.ambiguous_hits),
            ])


def write_hits_as_gff3(records: list[HitRecord], file: str) -> None:
    """
    Write a hit record as a GFF3 line, with the hit_class as an attribute.
    """
    with open(file, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for record in records:
            assert record.hit_class is not None, f"Hit class not assigned for hit {record.hit_id}"
            interval = record.hit
            writer.writerow([
                interval.seqid,
                interval.source,
                interval._type,
                str(interval.start),
                str(interval.end),
                interval.score,
                interval.strand,
                interval.phase,
                ";".join(f"{k}={v}" for k, v in interval.attrs.items())
            ])


def write_classified_gff3s(hits: dict[str, HitRecord], out_dir: str) -> None:
    """Partition hits by assigned class and write one GFF3 per class."""
    by_class: dict[HitClass, list[HitRecord]] = {hc: [] for hc in HitClass}
    for record in hits.values():
        assert record.hit_class is not None, f"unclassified hit {record.hit_id}"
        by_class[record.hit_class].append(record)
    for hit_class, filename in CLASS_TO_GFF3.items():
        write_hits_as_gff3(by_class[hit_class], os.path.join(out_dir, filename))


def main(argv=None) -> None:
    p = argparse.ArgumentParser(
        description="Classify bedtools -wao hits (exonic/intronic/antisense/intergenic/ambiguous); emit per-class GFF3s + gene CSV."
        )
    p.add_argument(
            "--loci",
            type=str,
            required=True,
            help="bedtools -wao output: hits vs loci"
            )
    p.add_argument(
            "--segments",
            type=str,
            required=True,
            help="bedtools -wao output: hits vs segments"
            )
    p.add_argument(
            "-o",
            "--out-dir",
            default=os.getcwd()
            )
    args = p.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)

    hits = collect_hits(args.loci, args.segments)
    assign_hit_classes(hits)                      # must precede gene summaries
    gene_summaries = collect_gene_summaries(hits)

    write_classified_gff3s(hits, args.out_dir)
    gene_summary_as_csv(gene_summaries,
                        file=os.path.join(args.out_dir, GENE_CENTRIC_CSV))

    counts: dict[HitClass, int] = {}
    for r in hits.values():
        counts[r.hit_class] = counts.get(r.hit_class, 0) + 1
    for hc in HitClass:
        print(f"{hc.value}: {counts.get(hc, 0)}", file=sys.stderr)


if __name__ == "__main__":
    main()


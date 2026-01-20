#!/usr/bin/env python3

import argparse
import csv
from dataclasses import dataclass, field
from enum import Enum
import sys
from typing import Optional


class ContextClass(Enum):
    EXON = "EXON"
    INTRON = "INTRON"
    OPPOSITE_STRAND_GENE = "OPPOSITE_STRAND_GENE"
    INTERGENIC = "INTERGENIC"
    AMBIGUOUS = "AMBIGUOUS"
    INCOMPLETE_MAPPING = "INCOMPLETE_MAPPING"


class AnnotationClass(Enum):
    GENE = "GENE"
    EXON = "EXON"
    OPP_GENE = "OPP_GENE"


def label_to_ann_class(label: str, *, default: Optional[AnnotationClass] = None) -> AnnotationClass:
    if label == "genes":
        return AnnotationClass.GENE
    if label == "exons":
        return AnnotationClass.EXON
    if label in ("opp_genes", "opposite_genes", "oppGenes"):
        return AnnotationClass.OPP_GENE
    if default is not None:
        return default
    raise ValueError(f"Unknown label {label!r} (expected genes/exons/opp_genes)")


@dataclass
class GffRecord:
    seqid: str
    source: str
    _type: str
    start_1based: int
    end_1based: int
    score: float | None
    strand: str
    phase: str
    attributes: dict[str, str]


def parse_attributes(attr: str) -> dict[str, str]:
    attr = attr.strip()
    if not attr or attr == ".":
        return {}
    out: dict[str, str] = {}
    for part in attr.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = ""
    return out


def parse_gff9(fields: list[str]) -> GffRecord:
    if len(fields) != 9:
        raise ValueError(f"Expected 9 fields for GFF record, got {len(fields)}: {fields}")
    seqid, source, _type, start, end, score, strand, phase, attrs = fields
    return GffRecord(
        seqid=seqid,
        source=source,
        _type=_type,
        start_1based=int(start),
        end_1based=int(end),
        score=None if score == "." else float(score),
        strand=strand,
        phase=phase,
        attributes=parse_attributes(attrs),
    )


@dataclass
class IntervalAgg:
    interval: GffRecord

    # same-strand
    same_strand_gene_ids: set[str] = field(default_factory=set)
    same_strand_exon_ids: set[str] = field(default_factory=set)
    same_strand_exon_gene_ids: set[str] = field(default_factory=set) # get parent gene IDs of same-strand exons

    # opposite-strand
    opp_strand_gene_ids: set[str] = field(default_factory=set)

    # incomplete mapping
    incomplete_exon2gene: bool = False

    # bookkeeping
    reasons: list[str] = field(default_factory=list)  # optional: for debugging


def read_gene_to_exon_map(path: str) -> dict[str, str]:
    """
    Expect a CSV with: gene_id,exon_id
    Lines starting with # ignored. Empty ignored.
    """
    m: dict[str, str] = {}
    with open(path, "r", newline="") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            gene_id, exon_id = parts[0].strip(), parts[1].strip()
            if exon_id and gene_id:
                m[exon_id] = gene_id
    return m


def iter_rows(path: str):
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            yield row


def has_complete_exon2gene_map(exon2gene: dict[str, str], id: str) -> bool:
    """Determine if id is gene or exon. Then return True if the mapping is complete."""
    genes = set(exon2gene.values())
    if id in genes:
        return True
    if id in exon2gene:
        return True

    return False


def update_from_intersections(
    aggs: dict[str, IntervalAgg],
    seen_any_hits: set[str],
    path: str,
    exon2gene: dict[str, str],
    is_opposite: bool,
    default_ann_class_if_no_label: Optional[AnnotationClass] = None,
) -> None:
    for row in iter_rows(path):
        if len(row) not in (18, 19):
            raise ValueError(f"{path}: expected 18/19 columns per row, got {len(row)}: {row}")
        if len(row) == 18:
            interval = parse_gff9(row[0:9])
            label = "opposite_genes" if is_opposite else "genes"
            annotation = parse_gff9(row[9:18])
        else:
            interval = parse_gff9(row[0:9])
            label = row[9]
            annotation = parse_gff9(row[10:19])

        # When using bedtools intersect with -wa -wb, the name that appears
        # in the label column (10th) is provided with the -names option.
        # Thus, we rely on that to determine the annotation class.
        ann_class = label_to_ann_class(label, default=default_ann_class_if_no_label)

        iid = interval.attributes.get("ID")
        aid = annotation.attributes.get("ID")
        if not iid or not aid:
            raise ValueError(f"{path}: missing ID. interval={interval.attributes} annotation={annotation.attributes}")

        # This interval has at least one hit, thus not intergenic
        seen_any_hits.add(iid)

        agg = aggs.get(iid)
        # Create if missing
        if agg is None:
            agg = aggs[iid] = IntervalAgg(interval=interval)

        # As an example of this, a pseudogene may have exons. These exons should not map to any gene.
        # Thus, we check for completeness of the exon2gene map.
        if not has_complete_exon2gene_map(exon2gene, aid):
            agg.incomplete_exon2gene = True

        # We are talking about intervals intersecting with genes here.
        if ann_class == AnnotationClass.GENE:
            if is_opposite:
                agg.opp_strand_gene_ids.add(aid)
            else:
                agg.same_strand_gene_ids.add(aid)

        # We are talking about intervals intersecting with exons here.
        elif ann_class == AnnotationClass.EXON:
            if is_opposite:
                print(f"Warning: saw EXON in opposite-strand file for interval {iid} (unexpected)", file=sys.stderr)
            agg.same_strand_exon_ids.add(aid)
            g = exon2gene.get(aid)
            if g:
                agg.same_strand_exon_gene_ids.add(g)
            else:
                # if exon->gene missing,
                agg.reasons.append(f"missing exon2gene for exon {aid}")

        elif ann_class == AnnotationClass.OPP_GENE:
            agg.opp_strand_gene_ids.add(aid)

        else:
            raise ValueError(f"Unhandled annotation class: {ann_class}")


def decide_context(agg: IntervalAgg) -> ContextClass:
    same_gene_identity = set(agg.same_strand_gene_ids) | set(agg.same_strand_exon_gene_ids) | set(agg.opp_strand_gene_ids)

    # Rule: ambiguous if intersects different genes (including exons from other genes)
    if len(same_gene_identity) > 1:
        return ContextClass.AMBIGUOUS

    # Determine presence of hits
    has_same_strand_gene = bool(agg.same_strand_gene_ids) or bool(agg.same_strand_exon_gene_ids)
    has_same_strand_exon = bool(agg.same_strand_exon_ids)
    has_opp_strand_gene = bool(agg.opp_strand_gene_ids)

    # If an exon does not map to any gene, we trust the map and thus, the exon may map to a pseudo-gene.
    if has_same_strand_exon and not has_same_strand_gene:
        # print to stderr for debugging
        agg.incomplete_exon2gene = True

    # EXON: intersects a gene AND an exon (and they resolve to same gene by identity)
    if has_same_strand_gene and has_same_strand_exon:
        return ContextClass.EXON

    # INTRON: intersects a gene but no exon
    if has_same_strand_gene and (not has_same_strand_exon):
        return ContextClass.INTRON

    # Opposite strand only
    if (not has_same_strand_gene) and (not has_same_strand_exon) and has_opp_strand_gene:
        return ContextClass.OPPOSITE_STRAND_GENE

    if agg.incomplete_exon2gene:
        return ContextClass.INCOMPLETE_MAPPING

    raise ValueError(f"Cannot classify interval {agg.interval.attributes.get('ID', '')}: no matching classification found")


def emit_csv_header(w: csv.writer) -> None:
    w.writerow(["interval_id", "context_class", "same_strand_gene_ids", "same_strand_exon_ids", "opp_strand_gene_ids"])


def emit_agg_row(w: csv.writer, iid: str, agg: IntervalAgg) -> None:
    same_genes = sorted(set(agg.same_strand_gene_ids) | set(agg.same_strand_exon_gene_ids))
    try:
        w.writerow([
            iid,
            decide_context(agg).value,
            ";".join(same_genes),
            ";".join(sorted(agg.same_strand_exon_ids)),
            ";".join(sorted(agg.opp_strand_gene_ids)),
        ])
    except ValueError as e:
        print(f"Warning: skipping interval {iid} due to error: {e}", file=sys.stderr)


def emit_nohits_as_intergenic(
    w: csv.writer,
    nohits_path: str,
    seen_any_hits: set[str],
) -> None:
    """
    Prints interval_id,INTERGENIC,,,, for intervals in nohits file
    that never appeared in any hit file.
    """
    for row in iter_rows(nohits_path):
        if len(row) != 9:
            raise ValueError(f"{nohits_path}: expected 9 columns, got {len(row)}: {row}")
        interval = parse_gff9(row)
        iid = interval.attributes.get("ID")
        if not iid:
            raise ValueError(f"{nohits_path}: interval missing ID attribute: {interval.attributes}")
        if iid in seen_any_hits:
            continue
        w.writerow([iid, ContextClass.INTERGENIC.value, "", "", ""])


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--same", required=True, help="Same-strand bedtools output (genes+exons)")
    ap.add_argument("--opp", required=False, help="Opposite-strand bedtools output (genes)")
    ap.add_argument("--nohits", required=True, help="GFF (9 cols) of intervals with no overlaps (bedtools -v output)")
    ap.add_argument("--gene2exon", required=True, help="CSV: gene_id,exon_id,")
    args = ap.parse_args()

    # Read gene->exon map, build exon->gene map
    gene2exon = read_gene_to_exon_map(args.gene2exon)

    # Datastructures definition
    aggs: dict[str, IntervalAgg] = {}
    seen_any_hits: set[str] = set()

    # Same-strand intersections
    update_from_intersections(
        aggs,
        seen_any_hits,
        args.same,
        gene2exon,
        is_opposite=False,
    )
    # Opposite-strand intersections
    update_from_intersections(
        aggs,
        seen_any_hits,
        args.opp,
        gene2exon,
        is_opposite=True,
        default_ann_class_if_no_label=AnnotationClass.OPP_GENE,
    )

    # Write output
    out_fh = sys.stdout
    w = csv.writer(out_fh, delimiter=",", lineterminator="\n")
    emit_csv_header(w)

    for iid in sorted(aggs.keys()):
        emit_agg_row(w, iid, aggs[iid])

    if args.nohits:
        emit_nohits_as_intergenic(w, args.nohits, seen_any_hits)


    return 0


if __name__ == "__main__":
    main()


include { AWK_MATCH_FEATURE_KEYWORD as AWK_MATCH_FEATURE_KEYWORD_loci       } from '../modules/local/awk/main'
include { AWK_MATCH_FEATURE_KEYWORD as AWK_MATCH_FEATURE_KEYWORD_segment    } from '../modules/local/awk/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_loci                     } from '../modules/local/bedtools/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_segment                  } from '../modules/local/bedtools/main'
include { PYTHON_JOIN_AND_CLASSIFY                                          } from '../modules/local/python/main'

workflow BUILD_LOCI_AND_SEGMENT_GFF{
    take:
    annotation  // tuple val(meta), path(annotation.gff3)

    main:

    // Gene features
    AWK_MATCH_FEATURE_KEYWORD_loci(
        annotation,
        params.kw_loci
    )

    // Exon features
    AWK_MATCH_FEATURE_KEYWORD_segment(
        annotation,
        params.kw_segment
    )

    emit:
    loci            = AWK_MATCH_FEATURE_KEYWORD_loci.out.matched_features
    segment         = AWK_MATCH_FEATURE_KEYWORD_segment.out.matched_features
}

workflow INTERSECT{
    take:
    tracks //val(meta), path(loci), path(segment) path(intervals)

    main:

    ch_loci = tracks.map { meta, loci, _segment, intervals ->
            [meta, intervals, loci]
        }
    ch_segment = tracks.map { meta, _loci, segment, intervals ->
            [meta, intervals, segment]
        }

    BEDTOOLS_INTERSECT_loci(
        ch_loci
    )

    BEDTOOLS_INTERSECT_segment(
        ch_segment
    )

    emit:
    loci_intersect      = BEDTOOLS_INTERSECT_loci.out.intersect_out
    segment_intersect   = BEDTOOLS_INTERSECT_segment.out.intersect_out
}

workflow JOIN_AND_CLASSIFY{
    take:
    join_in // val(meta), path(loci_intersect), path(segment_intersect)

    main:

    PYTHON_JOIN_AND_CLASSIFY(
        join_in
    )

    emit:
    gene_summary = PYTHON_JOIN_AND_CLASSIFY.out.gene_summary
    intervals    = PYTHON_JOIN_AND_CLASSIFY.out.intervals
}

include { BEDTOOLS_SORT as BEDTOOLS_SORT_gene                               } from '../../../modules/local/bedtools/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_transcript                         } from '../../../modules/local/bedtools/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_exon                               } from '../../../modules/local/bedtools/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_intervals                          } from '../../../modules/local/bedtools/main'

include { AWK_MATCH_FEATURE_KEYWORD as AWK_MATCH_FEATURE_KEYWORD_gene       } from '../../../modules/local/awk/main'
include { AWK_MATCH_FEATURE_KEYWORD as AWK_MATCH_FEATURE_KEYWORD_transcript } from '../../../modules/local/awk/main'
include { AWK_MATCH_FEATURE_KEYWORD as AWK_MATCH_FEATURE_KEYWORD_exon       } from '../../../modules/local/awk/main'

include { AWK_IDS_FROM_GFF                                                  } from '../../../modules/local/awk/main'


workflow NORMALIZE_AND_SORT_ANNOTATIONS {

    take:
    annotation // tuple val(meta), path(annotation.gff3)


    main:
    // Versions collector
    ch_versions = Channel.empty()

    // Gene features
    AWK_MATCH_FEATURE_KEYWORD_gene(
        annotation,
        params.kw_gene_feature,
        'genes'
    )
    ch_versions = ch_versions.mix(AWK_MATCH_FEATURE_KEYWORD_gene.out.versions)
    BEDTOOLS_SORT_gene(
        AWK_MATCH_FEATURE_KEYWORD_gene.out.matched_annotations,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT_gene.out.versions)

    // Transcript features
    AWK_MATCH_FEATURE_KEYWORD_transcript(
        annotation,
        params.kw_transcript_feature,
        'transcripts'
    )
    ch_versions = ch_versions.mix(AWK_MATCH_FEATURE_KEYWORD_transcript.out.versions)
    BEDTOOLS_SORT_transcript(
        AWK_MATCH_FEATURE_KEYWORD_transcript.out.matched_annotations,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT_transcript.out.versions)

    // Exon features
    AWK_MATCH_FEATURE_KEYWORD_exon(
        annotation,
        params.kw_exon_feature,
        'exons'
    )
    ch_versions = ch_versions.mix(AWK_MATCH_FEATURE_KEYWORD_exon.out.versions)
    BEDTOOLS_SORT_exon(
        AWK_MATCH_FEATURE_KEYWORD_exon.out.matched_annotations,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT_exon.out.versions)


    emit:
    genes       = BEDTOOLS_SORT_gene.out.sorted_gff
    transcripts = BEDTOOLS_SORT_transcript.out.sorted_gff
    exons       = BEDTOOLS_SORT_exon.out.sorted_gff
    versions    = ch_versions
}


workflow NORMALIZE_AND_SORT_INTERVALS {

    take:
    intervals  // tuple val(meta), path(intervals.gff3)


    main:
    // Versions collector
    ch_versions = Channel.empty()

    // Intervals sorting
    BEDTOOLS_SORT_intervals(
        intervals,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT_intervals.out.versions)
    AWK_IDS_FROM_GFF(
        BEDTOOLS_SORT_intervals.out.sorted_gff,
    )
    ch_versions = ch_versions.mix(AWK_IDS_FROM_GFF.out.versions)


    emit:
    intervals       = BEDTOOLS_SORT_intervals.out.sorted_gff
    interval_ids    = AWK_IDS_FROM_GFF.out.extracted_ids
    versions        = ch_versions
}

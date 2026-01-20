include { AWK_EXPORT_TRACKS_FROM_CLASSIFIED_INTERVALS   } from '../../../modules/local/awk/main'
include { AWK_MK_GENE_CENTRIC_TABLES                    } from '../../../modules/local/awk/main'
include { CONTEXT_FROM_BEDTOOLS_INTERSECT               } from '../../../modules/local/python/main'
include { COMM_BUILD_GENE_SETS                        } from '../../../modules/local/comm/main'


workflow CLASSIFY_CTX_CALLS {

    take:
    intervals_intersect             // tuple val(meta), path(intervals_intersect)
    intervals_opp_intersect         // tuple val(meta), path(intervals_opp_intersect)
    intervals_no_intersect          // tuple val(meta), path(intervals_no_intersect)
    gene_transcript_exon_hier_map   // tuple val(meta), path(gene_transcript_exon_hier_map.csv)
    intervals_sorted                // tuple val(meta), path(intervals_sorted)


    main:
    // Versions collector
    ch_versions = Channel.empty()

    CONTEXT_FROM_BEDTOOLS_INTERSECT(
        intervals_intersect,
        intervals_opp_intersect,
        intervals_no_intersect,
        gene_transcript_exon_hier_map
    )
    ch_versions = ch_versions.mix(CONTEXT_FROM_BEDTOOLS_INTERSECT.out.versions)

    AWK_EXPORT_TRACKS_FROM_CLASSIFIED_INTERVALS(
        CONTEXT_FROM_BEDTOOLS_INTERSECT.out.interval_context,
        intervals_sorted
    )
    ch_versions = ch_versions.mix(AWK_EXPORT_TRACKS_FROM_CLASSIFIED_INTERVALS.out.versions)

    AWK_MK_GENE_CENTRIC_TABLES(
        CONTEXT_FROM_BEDTOOLS_INTERSECT.out.interval_context,
    )
    ch_versions = ch_versions.mix(AWK_MK_GENE_CENTRIC_TABLES.out.versions)

    COMM_BUILD_GENE_SETS(
        gene_transcript_exon_hier_map,
        AWK_MK_GENE_CENTRIC_TABLES.out.gene_exonic_counts,
        AWK_MK_GENE_CENTRIC_TABLES.out.gene_intronic_counts,
        AWK_MK_GENE_CENTRIC_TABLES.out.gene_opp_strand_counts
    )
    ch_versions = ch_versions.mix(COMM_BUILD_GENE_SETS.out.versions)


    emit:
    versions = ch_versions
}

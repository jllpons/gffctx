include { BEDTOOLS_INTERSECT    } from '../../../modules/local/bedtools/main'


workflow INTERSECTIONS_CTX_CALLS {

    take:
    genes       // tuple val(meta), path(genes.sorted.gff3)
    exons       // tuple val(meta), path(exons.sorted.gff3)
    intervals   // tuple val(meta), path(intervals.sorted.gff3)


    main:
    // Versions collector
    ch_versions = Channel.empty()

    BEDTOOLS_INTERSECT(
        genes,
        exons,
        intervals,
        params.overlap_percent
        )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)


    emit:
    intervals_intersect     = BEDTOOLS_INTERSECT.out.intervals_intersect
    intervals_opp_intersect = BEDTOOLS_INTERSECT.out.intervals_opp_intersect
    intervals_no_intersect  = BEDTOOLS_INTERSECT.out.intervals_no_intersect
    versions                = ch_versions
}

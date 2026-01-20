include { NORMALIZE_AND_SORT_ANNOTATIONS        } from '../subworkflows/local/normalize_and_sort/main.nf'
include { NORMALIZE_AND_SORT_INTERVALS          } from '../subworkflows/local/normalize_and_sort/main.nf'
include { GENE_TRANSCRIPT_EXON_HIERARCHY_MAP    } from '../subworkflows/local/gene_transcript_exon_hier_map/main.nf'
include { INTERSECTIONS_CTX_CALLS               } from '../subworkflows/local/intersections_ctx_calls/main.nf'
include { CLASSIFY_CTX_CALLS                    } from '../subworkflows/local/classify_ctx_calls/main.nf'

include { generateQcMetadata                    } from '../modules/local/qc/main.nf'


workflow GFFCTX {


    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_annotation = samplesheet
        .map { meta, annotation, _intervals ->
            [meta, annotation]
        }

    ch_intervals  = samplesheet
        .map { meta, _annotation, intervals ->
            [meta, intervals]
        }

    // Step 0: Initialization and QC metadata generation
    generateQcMetadata()

    // Step 1: Normalize and sort GFF3 annotations & intervals
    NORMALIZE_AND_SORT_ANNOTATIONS(
        ch_annotation,
    )
    ch_versions = ch_versions.mix(NORMALIZE_AND_SORT_ANNOTATIONS.out.versions)
    NORMALIZE_AND_SORT_INTERVALS(
        ch_intervals,
    )
    ch_versions = ch_versions.mix(NORMALIZE_AND_SORT_INTERVALS.out.versions)

    // Step 2: Generate [gene -> transcript -> exon] hierarchy map from normalized & sorted GFF3
    GENE_TRANSCRIPT_EXON_HIERARCHY_MAP(
        NORMALIZE_AND_SORT_ANNOTATIONS.out.genes,
        NORMALIZE_AND_SORT_ANNOTATIONS.out.transcripts,
        NORMALIZE_AND_SORT_ANNOTATIONS.out.exons,
    )
    ch_versions = ch_versions.mix(GENE_TRANSCRIPT_EXON_HIERARCHY_MAP.out.versions)

    // Step 3: Perform intersections-based context calls
    INTERSECTIONS_CTX_CALLS(
        GENE_TRANSCRIPT_EXON_HIERARCHY_MAP.out.genes_filtered,
        GENE_TRANSCRIPT_EXON_HIERARCHY_MAP.out.exons_filtered,
        NORMALIZE_AND_SORT_INTERVALS.out.intervals,
    )
    ch_versions = ch_versions.mix(INTERSECTIONS_CTX_CALLS.out.versions)

    // Step 4: Give context to these intersections-based calls
    //         & classify them accordingly with ambiguity handling
    CLASSIFY_CTX_CALLS(
        INTERSECTIONS_CTX_CALLS.out.intervals_intersect,
        INTERSECTIONS_CTX_CALLS.out.intervals_opp_intersect,
        INTERSECTIONS_CTX_CALLS.out.intervals_no_intersect,
        GENE_TRANSCRIPT_EXON_HIERARCHY_MAP.out.gene_transcript_exon_hier_map,
        NORMALIZE_AND_SORT_INTERVALS.out.intervals,
    )
    ch_versions = ch_versions.mix(CLASSIFY_CTX_CALLS.out.versions)


}

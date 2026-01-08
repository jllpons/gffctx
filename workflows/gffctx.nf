include { NORMALIZE_AND_SORT                    } from '../subworkflows/local/normalize_and_sort/main.nf'
include { GENE_TRANSCRIPT_EXON_HIERARCHY_MAP    } from '../subworkflows/local/gene_transcript_exon_hier_map/main.nf'

include { generateQcMetadata                    } from '../modules/local/qc/main.nf'


workflow GFFCTX {


    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_annotation = samplesheet
        .map { meta, annotation  ->
            [meta, annotation]
        }
    /*
    ch_intervals  = samplesheet
        .map { meta, _annotation, intervals ->
            [meta, intervals]
        }
    */

    // Step 0: Initialization and QC metadata generation
    generateQcMetadata()

    // Step 1: Normalize and sort GFF3 annotations and intervals
    NORMALIZE_AND_SORT(
        ch_annotation,
    )
    ch_versions = ch_versions.mix(NORMALIZE_AND_SORT.out.versions)

    // Step 2: Generate [gene -> transcript -> exon] hierarchy map from normalized & sorted GFF3
    GENE_TRANSCRIPT_EXON_HIERARCHY_MAP(
        NORMALIZE_AND_SORT.out.genes,
        NORMALIZE_AND_SORT.out.transcripts,
        NORMALIZE_AND_SORT.out.exons,
    )
    ch_versions = ch_versions.mix(GENE_TRANSCRIPT_EXON_HIERARCHY_MAP.out.versions)


}

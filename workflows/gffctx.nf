include { BUILD_LOCI_AND_SEGMENT_GFF    } from '../subworkflows/main.nf'
include { INTERSECT                     } from '../subworkflows/main.nf'
include { PYTHON_JOIN_AND_CLASSIFY      } from '../subworkflows/main.nf'

workflow GFFCTX {

  take:
  samplesheet                                    // (sample, annotation, intervals)

  main:
    ch_versions = Channel.empty()

    ch_annotation = samplesheet.map { meta, annotation, _intervals -> [meta, annotation] }
    ch_intervals  = samplesheet.map { meta, _annotation, intervals  -> [meta, intervals] }

    BUILD_LOCI_AND_SEGMENT_GFF(ch_annotation)

    // Re-associate everything by meta so samples can't cross over
    ch_intersect_in = BUILD_LOCI_AND_SEGMENT_GFF.out.loci      // [meta, loci]
        .join(BUILD_LOCI_AND_SEGMENT_GFF.out.segment)          // [meta, loci, segment]
        .join(ch_intervals)                                    // [meta, loci, segment, intervals]

    INTERSECT(ch_intersect_in)

    // Re-associate everything by meta so samples can't cross over
    ch_join_in = INTERSECT.out.loci_intersect
        .join(INTERSECT.out.segment_intersect)                 // [meta, loci_intersect, segment_intersect]


    PYTHON_JOIN_AND_CLASSIFY(ch_join_in)
}

include { AWK_IDS_FROM_GFF                                                  } from '../../../modules/local/awk/main'
include { AWK_PARENTID_ID_FROM_GFF as AWK_PARENTID_ID_FROM_GFF_exons        } from '../../../modules/local/awk/main'
include { AWK_PARENTID_ID_FROM_GFF as AWK_PARENTID_ID_FROM_GFF_transcripts  } from '../../../modules/local/awk/main'
include { AWK_FILTER_GFF_BY_HIER_MAP_IDS                                    } from '../../../modules/local/awk/main'
include { BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON                   } from '../../../modules/local/python/main'


workflow GENE_TRANSCRIPT_EXON_HIERARCHY_MAP {

    take:
    genes       // tuple val(meta), path(genes.gff3)
    transcripts // tuple val(meta), path(transcripts.gff3)
    exons       // tuple val(meta), path(exons.gff3)


    main:
    // Versions collector
    ch_versions = Channel.empty()

    // Extract gene IDs
    AWK_IDS_FROM_GFF(genes)
    ch_versions = ch_versions.mix(AWK_IDS_FROM_GFF.out.versions)

    // Extract transcript IDs and their parent gene IDs
    AWK_PARENTID_ID_FROM_GFF_transcripts(transcripts)
    ch_versions = ch_versions.mix(AWK_PARENTID_ID_FROM_GFF_transcripts.out.versions)

    // Extract exon IDs and their parent transcript IDs
    AWK_PARENTID_ID_FROM_GFF_exons(exons)
    ch_versions = ch_versions.mix(AWK_PARENTID_ID_FROM_GFF_exons.out.versions)

    // Build hierarchical mapping gene -> transcript -> exon
    BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON(
        AWK_IDS_FROM_GFF.out.extracted_ids,
        AWK_PARENTID_ID_FROM_GFF_transcripts.out.extracted_parentid_ids,
        AWK_PARENTID_ID_FROM_GFF_exons.out.extracted_parentid_ids
    )
    ch_versions = ch_versions.mix(BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON.out.versions)


    AWK_FILTER_GFF_BY_HIER_MAP_IDS(
        genes,
        transcripts,
        exons,
        BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON.out.hier_map_gene_transcript_exon
    )
    ch_versions = ch_versions.mix(AWK_FILTER_GFF_BY_HIER_MAP_IDS.out.versions)

    emit:
    genes_filtered                  = AWK_FILTER_GFF_BY_HIER_MAP_IDS.out.filtered_genes
    transcripts_filtered            = AWK_FILTER_GFF_BY_HIER_MAP_IDS.out.filtered_transcripts
    exons_filtered                  = AWK_FILTER_GFF_BY_HIER_MAP_IDS.out.filtered_exons
    gene_transcript_exon_hier_map   = BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON.out.hier_map_gene_transcript_exon
    versions                        = ch_versions
}

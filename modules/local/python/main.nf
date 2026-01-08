process BUILD_HIERARCHICAL_MAPPING_GENE_TRANSCRIPT_EXON {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.12.12' }"

    publishDir "${params.outdir}/${sample_id.id}/intermediate", mode: 'copy', overwrite: true, pattern: "hier_map.gene-transcript-exon.csv"

    input:
    tuple val(sample_id), path(geneids)
    tuple val(sample_id_b), path(geneid_transcriptids)
    tuple val(sample_id_c), path(transcriptid_exonids)

    output:
    tuple val(sample_id), path("hier_map.gene-transcript-exon.csv"), emit: hier_map_gene_transcript_exon
    path 'versions.yml',                                             emit: versions

    script:
    """
    build_hierarchical_mapping_gene-transcript-exon.py --geneids ${geneids} \\
                                                        --geneid_transcriptids ${geneid_transcriptids} \\
                                                        --transcriptid_exonids ${transcriptid_exonids} \\
                                                        > hier_map.gene-transcript-exon.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}

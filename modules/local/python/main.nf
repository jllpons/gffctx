process PYTHON_JOIN_AND_CLASSIFY {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.11' :
    'quay.io/biocontainers/python:3.12.12' }"

    publishDir "${params.outdir}/${sample_id.id}/final", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bedtools_intersect_loci), path(bedtools_intersect_segments)

    output:
    path "gene_summary.csv",        emit: gene_summary
    path "*intervals.gff3",         emit: intervals
    //path 'versions.yml',            emit: versions

    script:
    """
    parse_bedtools_wao.py \\
        --loci ${bedtools_intersect_loci} \\
        --segment ${bedtools_intersect_segments} \\
    """
}


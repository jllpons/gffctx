process BEDTOOLS_INTERSECT {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    publishDir "${params.outdir}/${sample_id.id}/bedtools", mode: 'copy', overwrite: true, pattern: "*intersect.gff3"

    input:
    tuple val(sample_id), path(intervals), path(features)

    output:
    tuple val(sample_id), path("*.intersect.gff3"),     emit: intersect_out
    path 'versions.yml',                                emit: versions

    script:
    """
    intersect_output_file=${features.baseName}.intersect.gff3
    bedtools \\
        intersect \\
        -wao \\
        -a "$intervals" \\
        -b "$features" \\
        -f ${params.overlap_percent} \\
        >> "\$intersect_output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | awk '{print \$2}')
    END_VERSIONS
    """
}


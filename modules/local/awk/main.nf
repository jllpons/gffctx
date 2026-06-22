process AWK_MATCH_FEATURE_KEYWORD {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    publishDir "${params.outdir}/${sample_id.id}/annotation", mode: 'copy', overwrite: true, pattern: '*.gff3'

    input:
    tuple val(sample_id), path(annotation_file)
    val keyword

    output:
    tuple val(sample_id), path("${keyword}.gff3"),  emit: matched_features
    path 'versions.yml',                            emit: versions

    script:
    """
    output_file="${keyword}.gff3"
    timestamp=\$(date +"%Y-%m-%dT%H:%M:%S%z")

    echo "# gff_version: 3" > "\$output_file"
    echo "# md5_input_annotation: \$(md5sum "${annotation_file}" | awk '{print \$1}')" >> "\$output_file"
    echo "# feature_keyword_matching_cmd: gawk '!/^#/ && \\\$3 == \"${keyword}\"'" >> "\$output_file"
    echo "# operation_timestamp: \$timestamp" >> "\$output_file"
    gawk -F\$'\\t' '!/^#/ && \$3 == "${keyword}"' "${annotation_file}" >> "\$output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

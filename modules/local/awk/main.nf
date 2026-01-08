process AWK_MATCH_FEATURE_KEYWORD {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(sample_id), path(annotation_file)
    val keyword
    val filename_prefix

    output:
    tuple val(sample_id), path("${filename_prefix}.gff3"),  emit: matched_annotations
    path 'versions.yml',                                    emit: versions

    script:
    """
    output_file="${filename_prefix}.gff3"
    timestamp=\$(date +"%Y-%m-%dT%H:%M:%S%z")

    echo "# gff_version: 3" > "\$output_file"
    echo "# md5_input_annotation: \$(md5sum "${annotation_file}" | awk '{print \$1}')" >> "\$output_file"
    echo "# feature_keyword_matching_cmd: gawk '!/^#/ && \\\$3 == \"${keyword}\"'" >> "\$output_file"
    echo "# operation_timestamp: \$timestamp" >> "\$output_file"
    gawk '!/^#/ && \$3 == "${keyword}"' "${annotation_file}" >> "\$output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

process AWK_IDS_FROM_GFF {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(sample_id), path(gff_file)

    output:
    tuple val(sample_id), path("${gff_file.baseName}.IDs.csv"), emit: extracted_ids
    path 'versions.yml',                                        emit: versions

    script:
    """
    IDs_from_gff3.awk ${gff_file} | sort | uniq > "${gff_file.baseName}.IDs.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

process AWK_PARENTID_ID_FROM_GFF {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(sample_id), path(gff_file)

    output:
    tuple val(sample_id), path("${gff_file.baseName}.parentID_ID.csv"), emit: extracted_parentid_ids
    path 'versions.yml',                                                emit: versions

    script:
    """
    parentID_ID_from_gff3.awk ${gff_file} | sort -t',' -k1,1V -k2,2V | uniq > "${gff_file.baseName}.parentID_ID.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

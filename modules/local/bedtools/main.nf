process BEDTOOLS_SORT {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    publishDir "${params.outdir}/${sample_id.id}/intermediate", mode: 'copy', overwrite: true, pattern: "*.sorted.gff3"

    input:
    tuple val(sample_id), path(file)

    output:
    tuple val(sample_id), path("*.sorted.gff3"),    emit: sorted_gff
    path 'versions.yml',                            emit: versions

    script:
    """
    output_file="${file.baseName}.sorted.gff3"
    timestamp=\$(date +"%Y-%m-%dT%H:%M:%S%z")

    grep '^#' "$file" >> "\$output_file" || [ \$? -eq 1 ]
    echo "# md5_before_bedtools_sort: \$(md5sum "$file" | awk '{print \$1}')" >> "\$output_file"
    echo "# bedtools_sort_cmd: bedtools sort -i ${file}" >> "\$output_file"
    echo "# opperation_timestamp: \$timestamp" >> "\$output_file"

    bedtools sort -i "$file" >> "\$output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | awk '{print \$2}')
    END_VERSIONS
    """
}

process BEDTOOLS_INTERSECT {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    publishDir "${params.outdir}/${sample_id.id}/intermediate", mode: 'copy', overwrite: true, pattern: "*intersect.gff3"

    input:
    tuple val(sample_id),   path(genes)
    tuple val(sample_id_b), path(exons)
    tuple val(sample_id_c), path(intervals)
    val overlap_percent

    output:
    tuple val(sample_id), path("intervals.intersect.gff3"),     emit: intervals_intersect
    tuple val(sample_id), path("intervals.opp.intersect.gff3"), emit: intervals_opp_intersect
    tuple val(sample_id), path("intervals.no_intersect.gff3"),  emit: intervals_no_intersect
    path 'versions.yml',                                        emit: versions

    script:
    """
    intersect_output_file="intervals.intersect.gff3"
    opp_intersect_output_file="intervals.opp.intersect.gff3"
    no_intersect_output_file="intervals.no_intersect.gff3"
    timestamp=\$(date +"%Y-%m-%dT%H:%M:%S%z")

    echo "## Genes File Metadata ##" >> "\$intersect_output_file"
    grep '^#' "$genes" >> "\$intersect_output_file" || [ \$? -eq 1 ]
    echo "# md5_before_bedtools_intersect: \$(md5sum "$genes" | awk '{print \$1}')" >> "\$intersect_output_file"
    echo "## Exons File Metadata ##" >> "\$intersect_output_file"
    grep '^#' "$exons" >> "\$intersect_output_file" || [ \$? -eq 1 ]
    echo "# md5_before_bedtools_intersect: \$(md5sum "$exons" | awk '{print \$1}')" >> "\$intersect_output_file"
    echo "## Intervals File Metadata ##" >> "\$intersect_output_file"
    grep '^#' "$intervals" >> "\$intersect_output_file" || [ \$? -eq 1 ]
    echo "# md5_before_bedtools_intersect: \$(md5sum "$intervals" | awk '{print \$1}')" >> "\$intersect_output_file"
    echo "## Bedtools Intersect Operation Metadata  ##" >> "\$intersect_output_file"

    grep '^#' "\$intersect_output_file" >> "\$opp_intersect_output_file" || [ \$? -eq 1 ]
    grep '^#' "\$intersect_output_file" >> "\$no_intersect_output_file" || [ \$? -eq 1 ]


    echo "# bedtools_intersect_cmd: bedtools intersect -wa -wb -a ${intervals} -b ${genes} ${exons} -names genes exons -sorted -s -f ${overlap_percent}" >> "\$intersect_output_file"
    echo "# opperation_timestamp: \$timestamp" >> "\$intersect_output_file"
    bedtools \\
        intersect \\
        -wa \\
        -wb \\
        -a "$intervals" \\
        -b "$genes" "$exons" \\
        -names "genes" "exons" \\
        -sorted \\
        -s \\
        -f ${overlap_percent} \\
        >> "\$intersect_output_file"


    echo "# bedtools_intersect_cmd: bedtools intersect -wa -wb -a ${intervals} -b ${genes} -names opp_genes -sorted -S -f ${overlap_percent}" >> "\$opp_intersect_output_file"
    echo "# opperation_timestamp: \$timestamp" >> "\$opp_intersect_output_file"
    bedtools \\
        intersect \\
        -wa \\
        -wb \\
        -a "$intervals" \\
        -b "$genes" \\
        -names "opp_genes" \\
        -sorted \\
        -S \\
        -f ${overlap_percent} \\
        >> "\$opp_intersect_output_file"

    
    echo "# bedtools_intersect_cmd: bedtools intersect -a ${intervals} -b ${genes} -sorted -v -f ${overlap_percent}" >> "\$no_intersect_output_file"
    echo "# opperation_timestamp: \$timestamp" >> "\$no_intersect_output_file"
    bedtools \\
        intersect \\
        -a "$intervals" \\
        -b "$genes" \\
        -sorted \\
        -v \\
        -f ${overlap_percent} \\
        >> "\$no_intersect_output_file"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | awk '{print \$2}')
    END_VERSIONS
    """
}


process BEDTOOLS_SORT {
    tag "$sample_id.id"

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

    grep '^#' "$file" >> "\$output_file"
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


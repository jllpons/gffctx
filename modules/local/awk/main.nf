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
    IDs_from_gff3.awk ${gff_file} | sort -V | uniq > "${gff_file.baseName}.IDs.csv"

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

process AWK_EXPORT_TRACKS_FROM_CLASSIFIED_INTERVALS {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    publishDir "${params.outdir}/${sample_id.id}/tracks", mode: 'copy', overwrite: true, pattern: 'intervals.*.gff3'

    input:
    tuple val(sample_id),   path(context_file)
    tuple val(sample_id_b), path(tracks_file)

    output:
    tuple val(sample_id),   path('intervals.exonic.gff3'),              emit: exonic_intervals
    tuple val(sample_id_b), path('intervals.intronic.gff3'),            emit: intronic_intervals
    tuple val(sample_id),   path('intervals.opp_strand.gff3'),          emit: opp_strand_intervals
    tuple val(sample_id_b), path('intervals.intergenic.gff3'),          emit: intergenic_intervals
    tuple val(sample_id),   path('intervals.ambiguous.gff3'),           emit: ambiguous_intervals
    tuple val(sample_id_b), path('intervals.incomplete_mapping.gff3'),  emit: incomplete_mapping_intervals
    path 'versions.yml',    emit: versions

    script:
    """
    awk -F',' '/EXON/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.exonic.gff3 || true
    awk -F',' '/INTRON/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.intronic.gff3 || true
    awk -F',' '/OPPOSITE_STRAND_GENE/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.opp_strand.gff3 || true
    awk -F',' '/INTERGENIC/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.intergenic.gff3 || true
    awk -F',' '/AMBIGUOUS/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.ambiguous.gff3 || true
    awk -F',' '/INCOMPLETE_MAPPING/ {print "ID=" \$1 ";"}' ${context_file} | grep -Ff - ${tracks_file} > intervals.incomplete_mapping.gff3 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
        grep: \$(grep --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}


process AWK_FILTER_GFF_BY_HIER_MAP_IDS {
    tag "$sample_id.id"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
      ? 'docker://mbologna/docker-ripgrep:latest'
      : 'mbologna/docker-ripgrep:latest' }"

    publishDir "${params.outdir}/${sample_id.id}/intermediate", mode: 'copy', overwrite: true, pattern: 'filtered_*.gff3'

    input:
    tuple val(sample_id),   path(genes_file)
    tuple val(sample_id_b), path(transcripts_file)
    tuple val(sample_id_c), path(exons_file)
    tuple val(sample_id_d), path(hier_map)

    output:
    tuple val(sample_id),   path('filtered_genes.gff3'),        emit: filtered_genes
    tuple val(sample_id_b), path('filtered_transcripts.gff3'),  emit: filtered_transcripts
    tuple val(sample_id_c), path('filtered_exons.gff3'),        emit: filtered_exons
    path 'versions.yml',                                        emit: versions

    script:
    """
    awk -F',' '{print "ID=" \$1 ";"}' ${hier_map} | sort | uniq > hier_ids.genes.txt
    awk -F',' '{print "ID=" \$2 ";"}' ${hier_map} | sort | uniq > hier_ids.transcripts.txt
    awk -F',' '{print "ID=" \$3 ";"}' ${hier_map} | sort | uniq > hier_ids.exons.txt

    ( rg -Ff hier_ids.genes.txt ${genes_file} > filtered_genes.gff3 || true ) &
    ( rg -Ff hier_ids.transcripts.txt ${transcripts_file} > filtered_transcripts.gff3 || true ) &
    ( rg -Ff hier_ids.exons.txt ${exons_file} > filtered_exons.gff3 || true ) &
    wait

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
        rg: \$(rg -V | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
}


process AWK_MK_GENE_CENTRIC_TABLES {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    publishDir "${params.outdir}/${sample_id.id}", mode: 'copy', overwrite: true, pattern: "gene.*.csv"

    input:
    tuple val(sample_id), path(interval_context)

    output:
    tuple val(sample_id), path("gene.exonic.csv"),       emit: gene_exonic_counts
    tuple val(sample_id), path("gene.intronic.csv"),     emit: gene_intronic_counts
    tuple val(sample_id), path("gene.opp_strand.csv"),   emit: gene_opp_strand_counts
    path 'versions.yml',                                 emit: versions

    script:
    """
    grep 'EXON' ${interval_context} \
        | cut -d',' -f3 \
        | sort \
        | uniq -c \
        | gawk '{print \$2","\$1}' > gene.exonic.csv

    grep 'INTRON' ${interval_context} \
        | cut -d',' -f3 \
        | sort \
        | uniq -c \
        | gawk '{print \$2","\$1}' > gene.intronic.csv

    grep 'OPPOSITE_STRAND_GENE' ${interval_context} \
        | cut -d',' -f5 \
        | sort \
        | uniq -c \
        | gawk '{print \$2","\$1}' > gene.opp_strand.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

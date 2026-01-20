process COMM_BUILD_GENE_SETS {
    tag "$sample_id.id"

    publishDir "${params.outdir}/${sample_id.id}/gene_sets", mode: 'copy', overwrite: true, pattern: "*.list"

    input:
    tuple val(sample_id), path(hier_map_gene_transcript_exon)
    tuple val(sample_id_b), path(gene_exonic_counts)
    tuple val(sample_id_c), path(gene_intronic_counts)
    tuple val(sample_id_d), path(gene_opp_strand_counts)

    output:
    path 'genes_all.list',                                   emit: genes_all
    path 'genes_none.list',                                  emit: genes_none
    path 'genes_EIO.list',                                   emit: genes_EIO
    path 'genes_EI_only.list',                               emit: genes_EI_only
    path 'genes_EO_only.list',                               emit: genes_EO_only
    path 'genes_IO_only.list',                               emit: genes_IO_only
    path 'genes_E_only.list',                                emit: genes_E_only
    path 'genes_I_only.list',                                emit: genes_I_only
    path 'genes_O_only.list',                                emit: genes_O_only
    path 'versions.yml',                                    emit: versions

    script:
    """
    cat ${hier_map_gene_transcript_exon} | cut -d',' -f1 | sort | uniq > genes_all.list

    cat ${gene_exonic_counts} | cut -d',' -f1 | sort | uniq > genes_E.list
    cat ${gene_intronic_counts} | cut -d',' -f1 | sort | uniq > genes_I.list
    cat ${gene_opp_strand_counts} | cut -d',' -f1 | sort | uniq > genes_O.list

    cat genes_E.list genes_I.list genes_O.list | sort | uniq > genes_EIO_union.list
    comm -13 genes_EIO_union.list genes_all.list | sort | uniq > genes_none.list

    comm -12 genes_E.list genes_I.list | comm -12 - genes_O.list | sort | uniq > genes_EIO.list

    comm -12 genes_E.list genes_I.list | comm -23 - genes_O.list | sort | uniq > genes_EI_only.list
    comm -12 genes_E.list genes_O.list | comm -23 - genes_I.list | sort | uniq > genes_EO_only.list
    comm -12 genes_I.list genes_O.list | comm -23 - genes_E.list | sort | uniq > genes_IO_only.list

    comm -23 genes_EIO_union.list genes_I.list | comm -23 - genes_O.list | sort | uniq > genes_E_only.list
    comm -23 genes_EIO_union.list genes_E.list | comm -23 - genes_O.list | sort | uniq > genes_I_only.list
    comm -23 genes_EIO_union.list genes_E.list | comm -23 - genes_I.list | sort | uniq > genes_O_only.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comm: \$(comm --version 2>/dev/null | head -n1 | awk '{print \$4}' || echo "no version info")
        gawk: \$(gawk --version 2>/dev/null | awk 'NR==1{print \$3}' || echo "no version info")
    END_VERSIONS
    """
}

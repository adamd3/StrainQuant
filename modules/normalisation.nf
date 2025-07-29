process SUBSET_GENES {
    tag "$meta_merged"
    label 'process_medium'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/gene_counts", mode: 'symlink'

    input:
    path gpa_file
    path meta_merged
    val perc

    output:
    path 'gene_set_ST_*.tsv', emit: gene_subset

    script:
    def output_file = "gene_set_ST_${meta_merged.baseName}_${perc}perc.tsv"
    """
    subset_genes.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --perc=$perc \
        --ref_only=False \
        --outf=${output_file}
    """
}


process DESEQ_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/gene_counts", mode: 'symlink'

    input:
    tuple path(merged_counts), path(merged_lens)
    path gene_subset

    output:
    path 'deseq2_norm_counts.tsv', emit: norm_counts
    path 'deseq2_rpkm_counts.tsv', emit: rpkm_counts
    path 'deseq2_raw_counts.tsv', emit: scaled_counts

    script:
    def prefix = "deseq2"
    """
    DESeq2_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -p FALSE -t TRUE \
        -o ./ \
        --prefix ${prefix}
    """
}

process TMM_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/gene_counts", mode: 'symlink'

    input:
    tuple path(merged_counts), path(merged_lens)
    path gene_subset

    output:
    path 'tmm_norm_counts.tsv', emit: norm_counts
    path 'tmm_rpkm_counts.tsv', emit: rpkm_counts
    path 'tmm_raw_counts.tsv', emit: scaled_counts

    script:
    def prefix = "tmm"
    """
    TMM_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -p FALSE -t TRUE \
        -o ./ \
        --prefix ${prefix}
    """
}

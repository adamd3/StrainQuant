process SUBSET_GENES {
    meta {
        description = "Filters genes based on presence across strains to create a strain-typing gene set"
        keywords    = ["gene filtering", "strain typing", "gene presence", "subset"]
        authors     = ["@adamd3"]
        input       = [
            [ path(gpa_file), "Gene presence/absence file" ],
            [ path(meta_merged), "Merged sample metadata file" ],
            [ val(perc), "Minimum percentage threshold for gene presence across samples" ]
        ]
        output      = [
            [ path("gene_set_ST.tsv"), "Filtered gene set suitable for strain typing analysis" ]
        ]
    }

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
    meta {
        description = "Normalizes gene count data using DESeq2 size factors and calculates RPKM values"
        keywords    = ["normalization", "deseq2", "rpkm", "size factors", "gene expression"]
        authors     = ["@adamd3"]
        input       = [
            [ path(merged_counts), "Merged gene count matrix" ],
            [ path(merged_lens), "Gene length matrix" ],
            [ path(gene_subset), "Filtered gene set for analysis" ]
        ]
        output      = [
            [ path("norm_counts.tsv"), "DESeq2 normalized gene counts" ],
            [ path("rpkm_counts.tsv"), "RPKM normalized gene expression values" ],
            [ path("raw_counts.tsv"), "Raw gene counts scaled by library size" ]
        ]
    }

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
    meta {
        description = "Normalizes gene count data using Trimmed Mean of M-values (TMM) method and calculates RPKM values"
        keywords    = ["normalization", "tmm", "rpkm", "edger", "gene expression"]
        authors     = ["@adamd3"]
        input       = [
            [ path(merged_counts), "Merged gene count matrix" ],
            [ path(merged_lens), "Gene length matrix" ],
            [ path(gene_subset), "Filtered gene set for analysis" ]
        ]
        output      = [
            [ path("norm_counts.tsv"), "TMM normalized gene counts" ],
            [ path("rpkm_counts.tsv"), "RPKM normalized gene expression values" ],
            [ path("raw_counts.tsv"), "Raw gene counts scaled by library size" ]
        ]
    }

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

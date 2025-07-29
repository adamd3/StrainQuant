process UMAP_SAMPLES {
    tag "$norm_counts"
    label 'process_medium'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/umap_samples", mode: 'copy'

    input:
    path norm_counts
    path meta_merged
    val group

    output:
    path '*.{rds,png}', emit: umap_out

    script:
    """
    umap.R -n $norm_counts -m $meta_merged -g $group -o ./
    """
}

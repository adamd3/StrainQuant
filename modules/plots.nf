process UMAP_SAMPLES {
    meta {
        description = "Generates UMAP dimensionality reduction plots for sample visualization and clustering analysis"
        keywords    = ["umap", "dimensionality reduction", "visualization", "clustering", "samples"]
        authors     = ["@adamd3"]
        input       = [
            [ path(norm_counts), "Normalized gene count matrix" ],
            [ path(meta_merged), "Merged sample metadata file" ],
            [ val(group), "Grouping variable for sample coloring/annotation" ]
        ]
        output      = [
            [ path("*.{rds,png}"), "UMAP plots in PNG format and R data objects (RDS)" ]
        ]
    }

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

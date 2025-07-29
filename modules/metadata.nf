process CHECK_META_FILE {
    meta {
        description = "Validates and processes sample metadata file for strain-specific RNA-Seq analysis"
        keywords    = ["metadata", "validation", "preprocessing", "samples"]
        authors     = ["@adamd3"]
        input       = [
            [ path(metadata), "Sample metadata file in TSV format containing sample information" ]
        ]
        output      = [
            [ path("metadata_final.tsv"), "Validated and processed metadata file ready for downstream analysis" ]
        ]
    }

    tag "$metadata"
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata

    output:
    path 'metadata_final_*.tsv', emit: sample_metadata

    when:
    metadata.exists()

    script:
    def output_file = "metadata_final_${metadata.baseName}.tsv"
    """
    check_metadata.py $metadata ${output_file}
    """
}

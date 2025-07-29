process CHECK_META_FILE {
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

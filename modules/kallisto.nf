process KALLISTO_QUANT {
    meta {
        description = "Performs pseudo-alignment and quantification of RNA-Seq reads against strain-specific gene sequences using Kallisto"
        keywords    = ["quantification", "pseudo-alignment", "kallisto", "rna-seq", "strain-specific"]
        authors     = ["@adamd3"]
        input       = [
            [ path(gpa), "Gene presence/absence file defining strain-specific gene sets" ],
            [ val(meta), "Sample metadata map" ],
            [ path(reads), "Quality-trimmed FASTQ files" ],
            [ path(fasta), "Reference FASTA file containing gene sequences" ],
            [ val(strandedness), "Library strandedness (0=unstranded, 1=fr-stranded, 2=rf-stranded)" ]
        ]
        output      = [
            [ path("kallisto_*"), "Kallisto quantification output directory containing abundance estimates" ]
        ]
    }

    tag "$meta.sample_id"
    label 'process_high'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/kallisto_quant", mode: 'symlink'

    input:
    path gpa
    tuple val(meta), path(reads), path(fasta)
    val strandedness

    output:
    path "kallisto_${meta.sample_id}", emit: kallisto_out_dirs

    script:
    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (strandedness == 0){
        strand_arg = ""
    } else if (strandedness == 1){
        strand_arg = "--fr-stranded"
    } else if (strandedness==2){
        strand_arg = "--rf-stranded"
    }

    if (meta.paired_end) {
        // if trimming has not been performed, must symlink to match expected file names
        // kallisto params -l, -s are estimated from paired end data, but are required when using --single
        """
        [ ! -f  ${name}_1_val_1.fq.gz ] && ln -s ${reads[0]} ${name}_1_val_1.fq.gz
        [ ! -f  ${name}_2_val_2.fq.gz ] && ln -s ${reads[1]} ${name}_2_val_2.fq.gz
        kallisto index -i ${name}.kidx $fasta
        kallisto quant -t $task.cpus -i ${name}.kidx \
            ${strand_arg} -o kallisto_${name} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz
        """
    } else {
        """
        [ ! -f  ${name}_trimmed.fq.gz ] && ln -s $reads ${name}_trimmed.fq.gz
        kallisto index -i ${name}.kidx $fasta
        kallisto quant -t $task.cpus -i ${name}.kidx \
            ${strand_arg} --single -l ${params.fragment_len} -s ${params.fragment_sd} \
            -o kallisto_${name} ${name}_trimmed.fq.gz
        """
    }
}


process MERGE_COUNTS_AND_LENS {
    meta {
        description = "Merges Kallisto quantification results across all samples and extracts gene length information"
        keywords    = ["merge", "counts", "gene lengths", "kallisto", "matrix"]
        authors     = ["@adamd3"]
        input       = [
            [ path(gpa_file), "Gene presence/absence file" ],
            [ path(kallisto_dirs), "Collection of Kallisto output directories from all samples" ],
            [ path(meta_merged), "Merged sample metadata file" ]
        ]
        output      = [
            [ path("kallisto_merged_counts.tsv"), "Merged gene count matrix across all samples" ],
            [ path("kallisto_merged_lens.tsv"), "Gene length matrix for normalization" ]
        ]
    }

    tag "merge_counts_and_lens"
    label 'process_high'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/kallisto_quant", mode: 'symlink'

    input:
    path gpa_file
    path kallisto_dirs
    path meta_merged

    output:
    tuple path('kallisto_merged_counts_*.tsv'), path('kallisto_merged_lens_*.tsv'), emit: kallisto_merged_out

    script:
    def counts_file = "kallisto_merged_counts_${gpa_file.baseName}_${meta_merged.baseName}.tsv"
    def lens_file = "kallisto_merged_lens_${gpa_file.baseName}_${meta_merged.baseName}.tsv"
    """
    merge_kallisto_counts.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=${counts_file}

    merge_kallisto_lens.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=${lens_file}
    """
}

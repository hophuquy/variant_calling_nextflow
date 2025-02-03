#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TRIMMING_FASTP {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'staphb/fastp:latest' }"

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fastq_file1), path(fastq_file2)

    output:
    path "trimmed_${sample_name}_{1,2}.fastq.gz", emit: fastq
	path "trimming_fastp_report.html", emit: html

    publishDir "${params.outdir}/trimmed/data", mode: 'copy', overwrite: true

    script:
    """
    echo "Input file 1: $fastq_file1"
    echo "Input file 2: $fastq_file2"
    fastp -i $fastq_file1 -I $fastq_file2 --detect_adapter_for_pe \
        -o trimmed_${sample_name}_1.fastq.gz \
        -O trimmed_${sample_name}_2.fastq.gz \
        --html trimming_fastp_report.html
    """
}

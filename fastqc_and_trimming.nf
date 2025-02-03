#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fastq_file)

    output:
    path "*_fastqc.zip", emit: zip 
    path "*_fastqc.html", emit: html
	
    publishDir "${params.outdir}/untrimmed", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    # Run FastQC on forward and reverse read files
    fastqc $fastq_file
    """
}

process FASTQC_TRIMMED {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fastq_file)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
	
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    # Run FastQC on forward and reverse read files
    fastqc $fastq_file
    """
}

process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    tag "multiqc_report"

    input:
    path fastqc_files, stageAs: "fastq/*"

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
 	
    publishDir "${params.outdir}/untrimmed", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    echo "fastqc files: ${fastqc_files}"
    # Run MultiQC to aggregate FastQC results
    multiqc fastq/
    """
}

process MULTIQC_TRIMMED{
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    tag "multiqc_report"

    input:
    path fastqc_files, stageAs: "fastq/*"

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
 	
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    echo "fastqc files: ${fastqc_files}"
    # Run MultiQC to aggregate FastQC results
    multiqc fastq/
    """
}

process TRIMMING_FASTP {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'staphb/fastp:latest' }"

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fastq_file1), path(fastq_file2)

    output:
    path "${sample_name}_{1,2}.trimmed.fastq.gz", emit: fastq
	path "trimming_fastp_report.html", emit: html

    publishDir "${params.outdir}/trimmed/data", mode: 'copy', overwrite: true

    script:
    """
    echo "Input file 1: $fastq_file1"
    echo "Input file 2: $fastq_file2"
    fastp -i $fastq_file1 -I $fastq_file2 --detect_adapter_for_pe \
        -o ${sample_name}_1.trimmed.fastq.gz \
        -O ${sample_name}_2.trimmed.fastq.gz \
        --html trimming_fastp_report.html \
        -w ${task.cpus}
    """
}
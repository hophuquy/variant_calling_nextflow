#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process INDEX_REFERENCE_GENOME {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'staphb/samtools:latest' }"

    input:
    path ref

    output:
    path "${ref.baseName}.fa.fai" 

    script:
    """
    echo "STEP 1: Index Reference Genome"
    samtools faidx ${ref}
    """
}

process CREATE_REFERENCE_DICT {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path ref, stageAs: "ref/*"
    output:
    path "${ref[1].baseName}.dict" 
    script:
    """
    echo "STEP 2: Create Reference Dictionary"
    gatk CreateSequenceDictionary \
        -R ref/${ref[1].baseName}.fa \
        -O ${ref[1].baseName}.dict 
    """
}

process MARK_DUPLICATES  {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path sam_file

    output:
    path "${sam_file.baseName}.sorted_dedup_reads.bam" 

    script:
    """
    echo "STEP 3: Mark Duplicates and Sort - GATK4"
    gatk MarkDuplicatesSpark -I ${sam_file} -O ${sam_file.baseName}.sorted_dedup_reads.bam
    """
}

process BASE_QUALITY_RECALIBRATION {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path sorted_dedup_bam_tuple, stageAs: "data/*"

    output:
    path "${sorted_dedup_bam_tuple[0].baseName}_recal_data.table"

    publishDir "${params.outdir}/BaseRecalibrator_result", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory


    script:
    """
    echo "STEP 4: Base Quality Recalibration"
    gatk BaseRecalibrator \
        -I ${sorted_dedup_bam_tuple[0]} \
        -R ${sorted_dedup_bam_tuple[1].name} \
        --known-sites ${sorted_dedup_bam_tuple[4].name} \
        -O ${sorted_dedup_bam_tuple[0].baseName}_recal_data.table
    """
}

process APPLY_BQSR {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path sorted_dedup_bam_tuple, stageAs: "data/*"
    path recal_data_table 

    output:
    path "${sorted_dedup_bam_tuple[0].baseName}_sorted_dedup_bqsr_reads.bam" 

    script:
    """
    echo "STEP 4.2: Apply BQSR"
    gatk ApplyBQSR \
        -I ${sorted_dedup_bam_tuple[0]} \
        -R ${sorted_dedup_bam_tuple[1].name} \
        --bqsr-recal-file ${recal_data_table} \
        -O ${sorted_dedup_bam_tuple[0].baseName}_sorted_dedup_bqsr_reads.bam
    """
}

process COLLECT_METRICS {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path ref_tuple , stageAs: "data/*"
    path sorted_dedup_bqsr_bam_file

    output:
    path "${ref_tuple[0].baseName}_alignment_metrics.txt", emit: "align"
    path "${ref_tuple[0].baseName}_insert_size_metrics.txt", emit: "metrics"
    path "${ref_tuple[0].baseName}_insert_size_histogram.pdf", emit: "histogram"

    script:
    """
    echo "STEP 5: Collect Alignment & Insert Size Metrics"

    gatk CollectAlignmentSummaryMetrics \
        -R ${ref_tuple[1].name} \
        -I ${sorted_dedup_bqsr_bam_file} \
        -O ${ref_tuple[0].baseName}_alignment_metrics.txt

    gatk CollectInsertSizeMetrics \
        -I ${sorted_dedup_bqsr_bam_file} \
        --Histogram_FILE ${ref_tuple[0].baseName}_insert_size_histogram.pdf \
        -O ${ref_tuple[0].baseName}_insert_size_metrics.txt 
    """
}
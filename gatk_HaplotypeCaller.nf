#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CALL_VARIANTS_GERMLINE {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path ref_tuple , stageAs: "data/*"
    path sorted_dedup_bqsr_bam_file

    output:
    path "${ref_tuple[0].baseName}raw_variants.vcf"

    publishDir "${params.outdir}/gatk_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    echo "STEP 6: Call Variants - GATK HaplotypeCaller"
    gatk HaplotypeCaller --java-options "-Xmx64g -XX:ParallelGCThreads=${task.cpus}"  \
        -R ${ref_tuple[1].name} \
        -I ${sorted_dedup_bqsr_bam_file} \
        -O ${ref_tuple[0].baseName}raw_variants.vcf
        --native-pair-hmm-threads ${task.cpus}
    """
}

process EXTRACT_SNPS_AND_INDELS_GERMLINE {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'docker4sheng/gatk4:latest' }"

    input:
    path ref_tuple , stageAs: "data/*"
    path raw_variants

    output:
    path "${ref_tuple[0].baseName}_raw_snps.vcf"
    path "${ref_tuple[0].baseName}_raw_indels.vcf"

    publishDir "${params.outdir}/gatk_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    echo "STEP 6.1: Extract SNPs and INDELs"
    gatk SelectVariants \
        -R ${ref_tuple[1].name} \
        -V ${raw_variants} \
        --select-type SNP \
        -O ${ref_tuple[0].baseName}_raw_snps.vcf

    gatk SelectVariants \
        -R ${ref_tuple[1].name}  \
        -V ${raw_variants} \
        --select-type INDEL \
        -O ${ref_tuple[0].baseName}_raw_indels.vcf
    """
}
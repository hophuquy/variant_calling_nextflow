#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Process 1: Call Somatic Variants using Mutect2
process MUTECT2 {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"

    tag "Call somatic variants"
   

    input:
    path ref, stageAs: "data/*"
    path(tumor_bam)
    val(tumor_samplename) 
    path(normal_bam)
    val(normal_samplename) 


    output:
    path "${tumor_samplename}_somatic_variants_mutect2.vcf.gz*", emit: somatic_vcf
    path "${tumor_samplename}_f1r2.tar.gz", emit: f1r2_tar_gz


    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    script:
    """
    gatk Mutect2 -R ${ref[1]} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        --tumor ${tumor_samplename}.lite.1 \\
        -normal ${normal_samplename}.lite.1  \\
        --germline-resource ${ref[6]} \\
        --panel-of-normals ${ref[8]} \\
        -O ${tumor_samplename}_somatic_variants_mutect2.vcf.gz \\
        --f1r2-tar-gz ${tumor_samplename}_f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus}
    """
}

// Process 2: Estimate Cross-Sample Contamination (Tumor)
process GET_PILEUP_SUMMARIES {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"

    tag "Get pileup summaries for tumor"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    input:
    path ref, stageAs: "data/*"
    path bam
    path exon_region 
    val samplename

    output:
    path "${samplename}_getpileupsummaries.table"

    script:
    """
    mkdir -p tpm_workdir
    gatk GetPileupSummaries \\
        --java-options '-Xmx50G' --tmp-dir tpm_workdir \\
        -I ${bam} \\
        -V ${ref[6]} \\
        -L ${exon_region} \\
        -O ${samplename}_getpileupsummaries.table
    """
}

// Process 4: Calculate Contamination
process CALCULATE_CONTAMINATION {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"

    tag "Calculate contamination"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory


    input:
    path tumor_pileup_table
    path normal_pileup_table
    val tumor_samplename

    output:
    path "${tumor_samplename}_pair_calculatecontamination.table"

    script:
    """
    gatk CalculateContamination \\
        -I ${tumor_pileup_table} \\
        -matched ${normal_pileup_table} \\
        -O ${tumor_samplename}_pair_calculatecontamination.table
    """
}

// Process 5: Estimate Read Orientation Artifacts
process LEARN_READ_ORIENTATION_MODEL {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"

    tag "Learn read orientation model"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    input:
    path f1r2_tar_gz
    val samplename

    output:
    path "${samplename}_read-orientation-model.tar.gz"

    script:
    """
    gatk LearnReadOrientationModel \\
        -I ${f1r2_tar_gz} \\
        -O ${samplename}_read-orientation-model.tar.gz
    """
}

// Process 6: Filter Variants Called by Mutect2

process CREATE_FEATURE_INDEX {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"
    tag "create somatic vcf file index"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    input:
    path somatic_vcf

    output:
    path "${somatic_vcf.name}.tbi"

    script:
    """
    gatk IndexFeatureFile -I ${somatic_vcf} -O ${somatic_vcf.name}.tbi
    """
}

process FILTER_MUTECT_CALLS {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }"

    tag "Filter Mutect2 calls"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    input:
    path somatic_vcf, stageAs: "vcf/*" 
    path ref, stageAs: "data/*"
    path contamination_table
    path read_orientation_model
    val tumor_samplename
    output:
    path "${tumor_samplename}_somatic_variants_filtered_mutect2.vcf*"

    script:
    """
    

    gatk FilterMutectCalls \\
        -V ${somatic_vcf[0]} \\
        -R ${ref[1]} \\
        --contamination-table ${contamination_table} \\
        --ob-priors ${read_orientation_model} \\
        -O ${tumor_samplename}_somatic_variants_filtered_mutect2.vcf
    """
}

// Process 7: Annotate Variants using Funcotator
process FUNCOTATOR {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'broadinstitute/gatk:latest' }" 

    tag "Annotate variants with Funcotator"
    publishDir "${params.outdir}/gatk_mutect2_results", mode: 'copy', overwrite: true // Publish results to the 'fastqc_results' directory

    input:
    path filtered_vcf, stageAs: "vcf/*" 
    path ref, stageAs: "data/*"
    path data_sources
    val tumor_samplename

    output:
    path "${tumor_samplename}_somatic_variants_functotated.vcf"

    script:
    """
    gatk Funcotator \\
        --variant ${filtered_vcf[0]} \\
        --reference ${ref[1]} \\
        --ref-version hg38 \\
        --data-sources-path ${data_sources} \\
        --output ${tumor_samplename}_somatic_variants_functotated.vcf \\
        --output-file-format VCF
    """
}
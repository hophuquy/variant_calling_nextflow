#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BUILD_INDEX_BOWTIE2 {

    conda "bioconda::bowtie2=2.4.4 bioconda::samtools=1.16.1 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' :
        'biocontainers/bowtie2:v2.4.1_cv1' }"
    
    input:
    path reference_file // Path to the reference file

    output:
    path "index/*" // Path to store the Bowtie2 index files

    script:
    """
    mkdir index
    bowtie2-build ${reference_file} index/reference --threads ${task.cpus}
    """

}

process MAPPING_BOWTIE2 {
    conda "bioconda::bowtie2=2.4.4 bioconda::samtools=1.16.1 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' :
        'biocontainers/bowtie2:v2.4.1_cv1' }"

    input:
    path reference_index,  stageAs: "index/*"  // Staging the reference index as "index/reference"
    tuple val(sample_name), path(fastq_files)  // Input paired FASTQ files
    
  
    output:
    path "mapped/${sample_name}.sam" // Output SAM file for the mapped reads

    script:
    """
    #!/bin/bash
    mkdir -p mapped  // Ensure the output directory exists
    bowtie2 --verbose -x index/reference \
    -1 ${fastq_files[0]} \
    -2 ${fastq_files[1]} \
    --rg-id ${sample_name}.lite.1 \
    --rg "PL:ILLUMINA" \
    --rg "SM:${sample_name}.lite.1" \
    -S mapped/${sample_name}.sam \
    --threads 36
    """
}

process BUILD_INDEX_BWA {

     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0' :
        'biocontainers/bwa:v0.7.17_cv1' }"
    
    input:
    path reference_file // Path to the reference file

    output:
    path "*" // Path to store the Bowtie2 index files

    script:
    """
    mkdir index
    bwa index ${reference_file} 
    """
}

process MAPPING_BWA {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' :
        'biocontainers/bwa:v0.7.17_cv1' }"

    input:
    path reference_index,  stageAs: "index/*"  // Staging the reference index as "index/reference"
    tuple val(sample_name), path(fastq_files)  // Input paired FASTQ files
    
  
    output:
    path "mapped/${sample_name}.sam" // Output SAM file for the mapped reads

    script:
    """
    #!/bin/bash
    mkdir -p mapped  // Ensure the output directory exists
    
    # Run BWA MEM with proper escaping and line continuation
    bwa mem -t 36 -R "@RG\\tID:${sample_name}.lite.1\\tPL:ILLUMINA\\tSM:${sample_name}.lite.1" \\
    index/hg38.fa \\
    ${fastq_files[0]} \\
    ${fastq_files[1]} > mapped/${sample_name}.sam
    """
}


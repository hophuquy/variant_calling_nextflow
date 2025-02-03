#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC } from './fastqc_and_trimming.nf'
include { FASTQC_TRIMMED } from './fastqc_and_trimming.nf'
include { MULTIQC } from './fastqc_and_trimming.nf'
include { MULTIQC_TRIMMED} from './fastqc_and_trimming.nf'
include { TRIMMING_FASTP } from './fastqc_and_trimming.nf'
include { BUILD_INDEX_BWA } from './mapping.nf'
include { MAPPING_BWA } from './mapping.nf'
include { INDEX_REFERENCE_GENOME } from './gatk_preprocessing.nf'
include { CREATE_REFERENCE_DICT } from './gatk_preprocessing.nf'
include { MARK_DUPLICATES } from './gatk_preprocessing.nf'
include { BASE_QUALITY_RECALIBRATION } from './gatk_preprocessing.nf'
include { APPLY_BQSR } from './gatk_preprocessing.nf'
include { COLLECT_METRICS } from './gatk_preprocessing.nf' 
include { CALL_VARIANTS_GERMLINE  } from './gatk_HaplotypeCaller.nf'
include { EXTRACT_SNPS_AND_INDELS_GERMLINE } from './gatk_HaplotypeCaller.nf'
include { MUTECT2 } from './gatk_mutect2.nf'
include { GET_PILEUP_SUMMARIES } from './gatk_mutect2.nf'
include { CALCULATE_CONTAMINATION } from './gatk_mutect2.nf'
include {  LEARN_READ_ORIENTATION_MODEL  } from './gatk_mutect2.nf'
include {  CREATE_FEATURE_INDEX  } from './gatk_mutect2.nf'
include {  FILTER_MUTECT_CALLS  } from './gatk_mutect2.nf'
include {  FUNCOTATOR  } from './gatk_mutect2.nf'


Channel.fromPath(params.reads)
    .map { file -> tuple(file.simpleName, file) }
    .set { fastq_files }

Channel.fromFilePairs(params.reads, filePattern: '*_{1,2}.fastq', flat: true)
    .set { fastq_pair_files }

workflow {	
  
    // Quality control and trimming
    fastqc_results = FASTQC(fastq_files)
    multiqc_inputs = fastqc_results.zip.collect()

    // Run MultiQC after all FastQC processes are complete
    multiqc_results = MULTIQC(multiqc_inputs)

    // trimming
    trimmed_fastp_results = TRIMMING_FASTP(fastq_pair_files)

    trimmed_fastp_results.fastq.collect().flatten()
    .map { file -> 
        def fileName = file.simpleName.toString() // Ensure the file name is a string
        def sampleName = fileName.replaceAll(/_[12]/, '') // Remove _1.fastq.gz or _2.fastq.gz
        tuple(sampleName, file) // Create a tuple (sampleName, file)
    }
    .groupTuple() // Group by sampleName into a list of files
    .map { sampleName, files -> 
        tuple(sampleName, files.flatten()) // Sort the files to ensure R1 and R2 are in the correct order
    }
    .set { trimmed_fastq_pair_files } // Save the grouped pairs channel


    // trimmed_fastq_pair_files.view()

    // Create a channel from the list
    trimmed_fastps = trimmed_fastp_results.fastq.collect()
        // Wrap each element into a tuple
        .flatten()
        .map(file -> tuple(file.simpleName, file))

    // Quality control and trimming
    trimmed_fastqc_results = FASTQC_TRIMMED(trimmed_fastps)
    trimmed_multiqc_inputs = fastqc_results.zip.collect()

    // Run MultiQC after all FastQC processes are complete
    trimmed_multiqc_results = MULTIQC_TRIMMED(trimmed_multiqc_inputs)
    // Mapping 
    index_files = BUILD_INDEX_BWA(params.transcript)
    trimmed_fastq_pair_files.map{samplename, file -> tuple(samplename)}.collect().map{ it.sort({ a, b -> a <=> b }) }
    mapped_results = MAPPING_BWA(index_files, trimmed_fastq_pair_files)

    // GATK variant calling

    // Step 1: Index reference genome
    ref_fai = INDEX_REFERENCE_GENOME(params.transcript)
    // Ensure ref_fai is resolved to a file path
    ref_fai_resolved = ref_fai.map { fai -> tuple(fai)}
    ref_tuple = ref_fai.map { fai -> tuple(fai, params.transcript) }

    //  Step 2: Create reference dictionary
    ref_dict = CREATE_REFERENCE_DICT(ref_tuple)
    // Ensure ref_dict is resolved to a file path
    ref_dict_resolved = ref_fai.map { dict -> tuple(dict)}

    // Step 3 : mark deduplicates and Sort
    mark_deduplicates_result = MARK_DUPLICATES(mapped_results)
    // Combine channels into the tuple
    mark_deduplicates_result_tuple = mark_deduplicates_result.collect().map { it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten()
        .combine(ref_fai).combine(ref_dict)
        .map { file, fai, dict -> 
            tuple(
                file, 
                params.transcript, 
                fai,          // Resolved ref_fai
                dict,         // Resolved ref_dict
                params.known_sites, 
                params.known_sites_index,
                params.germline_resource,
                params.germline_resource_index,
                params.panel_of_normals,
                params.panel_of_normals_index
            )
        }
    // mark_deduplicates_result_tuple.view()
    // STEP 4: Base quality recalibration
    // # 1. build the model
    recalibration_results = BASE_QUALITY_RECALIBRATION(mark_deduplicates_result_tuple)
    recalibration_results_sorted = recalibration_results.collect().ifEmpty([]).map { it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten()

    // # 2. Apply the model to adjust the base quality scores
    apply_bqsr_result = APPLY_BQSR(mark_deduplicates_result_tuple, recalibration_results_sorted)
    apply_bqsr_files =apply_bqsr_result.collect().ifEmpty([]).map { it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten()

    // STEP 5: Collect Alignment & Insert Size Metrics
    metric_results = COLLECT_METRICS(mark_deduplicates_result_tuple, apply_bqsr_files)
    // apply_bqsr_files.collect().view()
    // // STEP 6: Call Variants - gatk haplotype caller
    // variant_calling_result = CALL_VARIANTS_GERMLINE(mark_deduplicates_result_tuple, apply_bqsr_files)
    // variant_calling_files = variant_calling_result.collect().ifEmpty([]).map { it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten()

    // // # extract SNPs & INDELS
    // EXTRACT_SNPS_AND_INDELS_GERMLINE(mark_deduplicates_result_tuple, variant_calling_files )

    // STEP 6: Call Variants - gatk mutect2

    normal_samplename = trimmed_fastq_pair_files.map{samplename, file -> tuple(samplename)}.collect().map{ it.sort({ a, b -> a <=> b }) }.flatten().filter( ~/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001/ )
    // Extract and filter tumor sample names
    tumor_samplename = trimmed_fastq_pair_files.map{samplename, file -> tuple(samplename)}.collect().map{ it.sort({ a, b -> a <=> b }) }.flatten().filter( ~/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001/ )
    normal_bam = apply_bqsr_files.filter( ~/.*HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.*/ )
    tumor_bam = apply_bqsr_files.filter( ~/.*HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.*/ )
    
    mutect2_result = MUTECT2(mark_deduplicates_result_tuple, tumor_bam, tumor_samplename, normal_bam, normal_samplename)

    // STEP 7: Estimate cross-sample contamination
    sample_names = trimmed_fastq_pair_files.map{samplename, file -> tuple(samplename)}.collect().map{ it.sort({ a, b -> a <=> b }) }.flatten()
    apply_bqsr_files.view()
    sample_names.view()
    pileup_summaries_results = GET_PILEUP_SUMMARIES(mark_deduplicates_result_tuple, apply_bqsr_files, params.exon_region, sample_names)
    tumor_pileup_summaries_files = pileup_summaries_results.collect().map{ it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten().filter( ~/.*HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.*/ )
    normal_pileup_summaries_files = pileup_summaries_results.collect().map{ it.sort({ a, b -> a.simpleName <=> b.simpleName }) }.flatten().filter( ~/.*HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.*/ )

    // Calculate contamination
    contamination_result = CALCULATE_CONTAMINATION(tumor_pileup_summaries_files, normal_pileup_summaries_files, tumor_samplename)
    contamination_result.view()

    // STEP 8: Estimate read orientation artifacts
    // read orientation
    orientation_model_result = LEARN_READ_ORIENTATION_MODEL(mutect2_result.f1r2_tar_gz, tumor_samplename)
    orientation_model_result.view()

    // STEP 9: Filter Variants Called By Mutect2
    somatic_mutation_filter_result = FILTER_MUTECT_CALLS(mutect2_result.somatic_vcf, mark_deduplicates_result_tuple, contamination_result, orientation_model_result, tumor_samplename)
    somatic_mutation_filter_result.view()
    funcotator_result = FUNCOTATOR(somatic_mutation_filter_result, mark_deduplicates_result_tuple, params.annotaion_database_dir, tumor_samplename)
}


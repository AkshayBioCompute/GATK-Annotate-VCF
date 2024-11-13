#!/usr/bin/env nextflow

// Define paths
def ref = "/home/akshay/Akshay/WGS/hg38.fa.gz"
def results = "/home/akshay/Akshay/WGS/results"
def raw_snps_vcf = "/home/akshay/Akshay/WGS/results/raw_snps.vcf"
def raw_indels_vcf = "/home/akshay/Akshay/WGS/results/raw_indels.vcf"
def funcotator_data_sources = "/home/akshay/Akshay/tools/funcotator/hg38/funcotator_dataSources.v1.7.20200521g"

process filter_variants {
    input:
    path raw_snps_vcf
    path raw_indels_vcf
    path ref

    output:
    path "${results}/filtered_snps.vcf"
    path "${results}/filtered_indels.vcf"

    script:
    """
    gatk VariantFiltration \
        -R ${ref} \
        -V ${raw_snps_vcf} \
        -O ${results}/filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"

    gatk VariantFiltration \
        -R ${ref} \
        -V ${raw_indels_vcf} \
        -O ${results}/filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"
    """
}

process select_variants {
    input:
    path filtered_snps_vcf
    path filtered_indels_vcf

    output:
    path "${results}/analysis-ready-snps.vcf"
    path "${results}/analysis-ready-indels.vcf"

    script:
    """
    gatk SelectVariants \
        --exclude-filtered \
        -V ${filtered_snps_vcf} \
        -O ${results}/analysis-ready-snps.vcf

    gatk SelectVariants \
        --exclude-filtered \
        -V ${filtered_indels_vcf} \
        -O ${results}/analysis-ready-indels.vcf
    """
}

process annotate_variants {
    input:
    path snps_vcf
    path indels_vcf
    path ref
    path funcotator_data_sources

    output:
    path "${results}/analysis-ready-snps-functotated.vcf"
    path "${results}/analysis-ready-indels-functotated.vcf"

    script:
    """
    gatk Funcotator \
        --variant ${snps_vcf} \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path ${funcotator_data_sources} \
        --output ${results}/analysis-ready-snps-functotated.vcf \
        --output-file-format VCF

    gatk Funcotator \
        --variant ${indels_vcf} \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path ${funcotator_data_sources} \
        --output ${results}/analysis-ready-indels-functotated.vcf \
        --output-file-format VCF
    """
}

workflow {
    filter_variants()
    select_variants()
    annotate_variants()
}


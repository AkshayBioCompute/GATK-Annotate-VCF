import os

# Directories
ref = "/home/akshay/Akshay/WGS/hg38.fa.gz"
results = "/home/akshay/Akshay/WGS/results"
raw_snps_vcf = "/home/akshay/Akshay/WGS/results/raw_snps.vcf"
raw_indels_vcf = "/home/akshay/Akshay/WGS/results/raw_indels.vcf"
funcotator_data_sources = "/home/akshay/Akshay/tools/funcotator/hg38/funcotator_dataSources.v1.7.20200521g"

# Rule to filter variants
rule filter_variants:
    input:
        snps=raw_snps_vcf,
        indels=raw_indels_vcf,
        ref=ref
    output:
        filtered_snps = os.path.join(results, "filtered_snps.vcf"),
        filtered_indels = os.path.join(results, "filtered_indels.vcf")
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.snps} \
            -O {output.filtered_snps} \
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
            -R {input.ref} \
            -V {input.indels} \
            -O {output.filtered_indels} \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 200.0" \
            -filter-name "SOR_filter" -filter "SOR > 10.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"
        """

# Rule to select PASS variants
rule select_variants:
    input:
        filtered_snps = os.path.join(results, "filtered_snps.vcf"),
        filtered_indels = os.path.join(results, "filtered_indels.vcf")
    output:
        snps_pass = os.path.join(results, "analysis-ready-snps.vcf"),
        indels_pass = os.path.join(results, "analysis-ready-indels.vcf")
    shell:
        """
        gatk SelectVariants \
            --exclude-filtered \
            -V {input.filtered_snps} \
            -O {output.snps_pass}

        gatk SelectVariants \
            --exclude-filtered \
            -V {input.filtered_indels} \
            -O {output.indels_pass}
        """

# Rule to annotate variants with Funcotator
rule annotate_variants:
    input:
        snps_pass = os.path.join(results, "analysis-ready-snps.vcf"),
        indels_pass = os.path.join(results, "analysis-ready-indels.vcf"),
        ref = ref,
        funcotator_data_sources = funcotator_data_sources
    output:
        snps_annotated = os.path.join(results, "analysis-ready-snps-functotated.vcf"),
        indels_annotated = os.path.join(results, "analysis-ready-indels-functotated.vcf")
    shell:
        """
        gatk Funcotator \
            --variant {input.snps_pass} \
            --reference {input.ref} \
            --ref-version hg38 \
            --data-sources-path {input.funcotator_data_sources} \
            --output {output.snps_annotated} \
            --output-file-format VCF

        gatk Funcotator \
            --variant {input.indels_pass} \
            --reference {input.ref} \
            --ref-version hg38 \
            --data-sources-path {input.funcotator_data_sources} \
            --output {output.indels_annotated} \
            --output-file-format VCF
        """


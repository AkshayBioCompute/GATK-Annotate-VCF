# File paths
REFERENCE_GENOME = "/home/akshay/Akshay/WGS/hg38.fa.gz"
RAW_SNPS_VCF = "/home/akshay/Akshay/WGS/results/raw_snps.vcf"
RAW_INDELS_VCF = "/home/akshay/Akshay/WGS/results/raw_indels.vcf"
FUNCOTATOR_DATA_SOURCES = "/home/akshay/Akshay/tools/funcotator/hg38/funcotator_dataSources.v1.7.20200521g"

# Filter criteria (e.g., for VariantFiltration)
FILTER_CRITERIA = {
    "QD": 2.0,
    "FS": 60.0
}

# Output directory
OUTPUT_DIR = "/home/akshay/Akshay/WGS/results"

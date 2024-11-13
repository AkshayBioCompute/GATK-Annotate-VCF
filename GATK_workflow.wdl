workflow VariantFilteringAnnotation {

  input {
    File ref                     # Reference genome file
    File raw_snps                # Raw SNPs VCF file
    File raw_indels              # Raw INDELs VCF file
    String results               # Directory for output files
    Map[String, Float] snp_filters
    Map[String, Float] indel_filters
    Map[String, Int] genotype_filters
    String funcotator_data_sources # Path to Funcotator data sources
    String ref_version           # Reference version for Funcotator (e.g., "hg38")
  }

  # Step 1: Filter SNP Variants
  call FilterVariants as FilterSNPs {
    input:
      ref = ref,
      input_vcf = raw_snps,
      output_vcf = "${results}/filtered_snps.vcf",
      filters = snp_filters,
      genotype_filters = genotype_filters
  }

  # Step 2: Filter INDEL Variants
  call FilterVariants as FilterINDELs {
    input:
      ref = ref,
      input_vcf = raw_indels,
      output_vcf = "${results}/filtered_indels.vcf",
      filters = indel_filters,
      genotype_filters = genotype_filters
  }

  # Step 3: Select Variants that PASS filters
  call SelectVariants {
    input:
      filtered_snps = FilterSNPs.output_vcf,
      filtered_indels = FilterINDELs.output_vcf,
      output_snps = "${results}/analysis-ready-snps.vcf",
      output_indels = "${results}/analysis-ready-indels.vcf"
  }

  # Step 4: Annotate Variants with Funcotator
  call Funcotator as AnnotateSNPs {
    input:
      variant = SelectVariants.output_snps,
      reference = ref,
      ref_version = ref_version,
      data_sources = funcotator_data_sources,
      output_vcf = "${results}/analysis-ready-snps-functotated.vcf"
  }
  
  call Funcotator as AnnotateINDELs {
    input:
      variant = SelectVariants.output_indels,
      reference = ref,
      ref_version = ref_version,
      data_sources = funcotator_data_sources,
      output_vcf = "${results}/analysis-ready-indels-functotated.vcf"
  }
}

task FilterVariants {
  input {
    File ref
    File input_vcf
    String output_vcf
    Map[String, Float] filters
    Map[String, Int] genotype_filters
  }

  command <<<
    gatk VariantFiltration \
      -R ~{ref} \
      -V ~{input_vcf} \
      -O ~{output_vcf} \
      $(for key, value in filters select("-filter-name", key + "_filter", "-filter", key + " < " + value))
      $(for key, value in genotype_filters select("-genotype-filter-expression", key + " < " + value, "-genotype-filter-name", key + "_filter"))
  >>>
  
  output {
    File output_vcf
  }
}

task SelectVariants {
  input {
    File filtered_snps
    File filtered_indels
    String output_snps
    String output_indels
  }

  command <<<
    gatk SelectVariants --exclude-filtered -V ~{filtered_snps} -O ~{output_snps}
    gatk SelectVariants --exclude-filtered -V ~{filtered_indels} -O ~{output_indels}
  >>>

  output {
    File output_snps
    File output_indels
  }
}

task Funcotator {
  input {
    File variant
    File reference
    String ref_version
    String data_sources
    String output_vcf
  }

  command <<<
    gatk Funcotator \
      --variant ~{variant} \
      --reference ~{reference} \
      --ref-version ~{ref_version} \
      --data-sources-path ~{data_sources} \
      --output ~{output_vcf} \
      --output-file-format VCF
  >>>

  output {
    File output_vcf
  }
}


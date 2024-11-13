import subprocess

def annotate_variants(vcf_file, ref_genome, funcotator_data_sources, output_file):
    """
    Annotates variants in a VCF file using GATK Funcotator.

    :param vcf_file: Path to the VCF file to annotate
    :param ref_genome: Path to the reference genome file
    :param funcotator_data_sources: Path to Funcotator data sources
    :param output_file: Path to save the annotated VCF
    """
    command = [
        "gatk", "Funcotator",
        "--variant", vcf_file,
        "--reference", ref_genome,
        "--ref-version", "hg38",
        "--data-sources-path", funcotator_data_sources,
        "--output", output_file,
        "--output-file-format", "VCF"
    ]
    subprocess.run(command, check=True)

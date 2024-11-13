import vcfpy

def filter_variants(vcf_file, filter_criteria):
    """
    Filters variants in a VCF file based on provided filter criteria.

    :param vcf_file: Path to the input VCF file
    :param filter_criteria: Dictionary of filter rules (e.g., QD < 2.0, FS > 60.0)
    :return: Filtered VCF object
    """
    reader = vcfpy.Reader.from_path(vcf_file)
    filtered_variants = []
    for record in reader:
        if passes_filters(record, filter_criteria):
            filtered_variants.append(record)
    return filtered_variants

def passes_filters(record, criteria):
    """
    Checks if a variant record passes the provided filter criteria.
    """
    for filter_name, filter_value in criteria.items():
        if filter_name == "QD" and record.QUAL < filter_value:
            return False
        if filter_name == "FS" and record.INFO["FS"] > filter_value:
            return False
    return True

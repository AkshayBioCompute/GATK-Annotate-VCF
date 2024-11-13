import os

def check_file_exists(file_path):
    """
    Checks if a given file exists.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    return True

def read_vcf(file_path):
    """
    Reads a VCF file and returns a list of records.
    """
    with open(file_path, 'r') as vcf_file:
        return vcf_file.readlines()

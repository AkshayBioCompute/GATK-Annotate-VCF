# GATK-Annotate-VCF

This repository contains a comprehensive pipeline for genomic data analysis, including tasks like quality control, read alignment, variant calling, annotation, and filtering. The pipeline supports multiple workflow engines, including Snakemake, Nextflow, and WDL. This flexibility allows the user to choose the most suitable execution environment for their project.

## Table of Contents

1. [Introduction](#introduction)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [File Structure](#file-structure)
5. [Usage](#usage)
6. [Execution with Snakemake](#execution-with-snakemake)
7. [Execution with Nextflow](#execution-with-nextflow)
8. [Execution with WDL](#execution-with-wdl)
9. [License](#license)

---

## Introduction

This pipeline is designed for analyzing high-throughput genomic sequencing data, specifically focusing on:

- **Quality Control** using FastQC
- **Mapping** sequencing reads to a reference genome using BWA
- **Marking Duplicates** and sorting the reads
- **Base Quality Recalibration** using GATK
- **Variant Calling and Filtering**
- **Variant Annotation** using GATK Funcotator

---

## Requirements

### Software Dependencies

- Python 3.7+
- Snakemake (for Snakemake workflow)
- Nextflow (for Nextflow workflow)
- Cromwell (for WDL workflow)
- GATK (Genome Analysis Toolkit)
- BWA (Burrows-Wheeler Aligner)
- FastQC
- GATK Funcotator

### Python Libraries

- pandas
- numpy
- scipy
- vcfpy
- scikit-allel
- tqdm
- pyyaml
- requests

To install the required Python libraries, you can use the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
```

### System Requirements

- Sufficient computational resources for genome analysis (e.g., 16GB RAM and multi-core processors recommended)
- Storage space for raw sequencing data, reference genomes, and intermediate files

---

## Installation

Clone this repository to your local machine:

```bash
git clone https://github.com/AkshayBioCompute/GATK-Annotate-VCF.git
cd GATK-Annotate-VCF.git
```

Install the required Python dependencies:

```bash
pip install -r requirements.txt
```

---

## File Structure

```plaintext
GATK-Annotate-VCF/
├── README.md                       # This readme file
├── requirements.txt                # Python dependencies
├── src/                            # Python scripts for the pipeline
│   ├── __init__.py
│   ├── data_processing.py
│   ├── variant_annotation.py
│   ├── utils.py
│   ├── pipeline_config.py
│   └── logger.py
├── snakefile                       # Snakemake workflow file
├── nextflow.config                 # Nextflow config file
├── main.nf                         # Nextflow workflow file
├── WDL/
│   └── pipeline.wdl                # WDL workflow file
├── input.json                      # JSON input file for WDL
└── results/                        # Output directory for analysis results
```

---

## Usage

### 1. Configure the Pipeline

Before running the pipeline, you must configure the paths and parameters in the configuration files.

- **For Snakemake:** Update the paths in the `Snakefile` and `src/pipeline_config.py`.
- **For Nextflow:** Update the paths in the `nextflow.config` file and `main.nf`.
- **For WDL:** Set the paths in the `input.json` file and `pipeline.wdl`.

### 2. Running the Pipeline

#### Execution with Snakemake

1. Create a virtual environment and install dependencies:

    ```bash
    conda create -n genome-pipeline snakemake
    conda activate genome-pipeline
    pip install -r requirements.txt
    ```

2. Run the pipeline:

    ```bash
    snakemake -s Snakefile --cores <number_of_cores>
    ```

#### Execution with Nextflow

1. Install Nextflow:

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

2. Run the pipeline:

    ```bash
    nextflow run main.nf -c nextflow.config
    ```

#### Execution with WDL

1. Install Cromwell (or use an existing instance):

    - Download Cromwell from [Cromwell GitHub Releases](https://github.com/broadinstitute/cromwell/releases).
    - Run the Cromwell server.

2. Submit the WDL workflow:

    ```bash
    java -jar cromwell.jar run WDL/pipeline.wdl -i input.json
    ```

---

## Execution Flow

The pipeline follows these main steps:

1. **Quality Control (QC):** Run FastQC to check the quality of the sequencing reads.
2. **Read Alignment:** Use BWA to align reads to the reference genome.
3. **Mark Duplicates and Sort:** Mark duplicates and sort the BAM file to prepare for variant calling.
4. **Base Quality Recalibration:** Perform base quality recalibration using GATK to improve accuracy in downstream analyses.
5. **Variant Calling and Filtering:** Use GATK HaplotypeCaller to call variants, followed by variant filtering.
6. **Variant Annotation:** Annotate the variants using GATK Funcotator.

Each step can be customized by modifying the corresponding configuration files.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

This `README.md` gives an overview of how to set up, configure, and run the genome analysis pipeline using different workflow engines. Make sure to adjust paths and configuration settings to match your local setup or computational environment.

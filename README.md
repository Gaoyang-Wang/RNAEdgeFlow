# RNAEdgeFlow

RNAEdgeFlow is a modular pipeline for processing RNA 5′/3′ end-tag and internal reads.  
It supports end-tag classification (5P, 3P, internal), UMI handling, adapter trimming, alignment, expression quantification, coverage profiling, and metagene analysis.

---

## Features

- Automatic separation of 5P, 3P, and internal reads from paired-end FASTQ.  
- UMI extraction and deduplication (via umi_tools).  
- Adapter trimming and quality control (via cutadapt and fastp).  
- Alignment with STAR and expression quantification with StringTie.  
- Coverage profiling and base content statistics.  
- Two major functions:  
  1. **pipeline** → full end-to-end pipeline for a sample, from FASTQ to processed BAMs, QC, expression matrices, and coverage.  
  2. **metagene** → metagene coverage plots (e.g., TSS/TES) from BAMs or directories of BAM files.  

---

## Installation

Step 1. Clone the repository  
git clone https://github.com/Gaoyang-Wang/RNAEdgeFlow.git  
cd RNAEdgeFlow  

Step 2. Create the Conda environment  
bash install.sh  

This will create a Conda environment named `rnaedgeflow` and install all required dependencies (Python, R, Bioconductor packages, bedtools, samtools, STAR, fastp, cutadapt, umi_tools, etc.).  

Step 3. Activate the environment  
conda activate rnaedgeflow  

---

## Usage

### 1. Run the full pipeline
./rnaedgeflow pipeline --profile path/to/sample.profile  

`--profile` points to a `.profile` file describing all required parameters (sample name, input FASTQ prefix, STAR genome index, barcode/UMI patterns, etc.).  

Outputs will be created under:  
OutputDir/SampleName/  
  ├── process/           # intermediate FASTQs and BAMs  
  ├── result/            # final outputs  
  │   ├── stat/          # QC and read statistics  
  │   ├── internal_expr/ # StringTie expression  
  │   └── terminal_bed/  # coverage profiles  
  └── log/               # pipeline logs  

Example:  
./rnaedgeflow pipeline --profile example_data/example.profile  

### 2. Metagene profiling

(a) Scan all BAMs in a directory  
./rnaedgeflow metagene --inputDir path/to/bam_dir --outputDir results/metagene_out  

(b) Explicitly provide BAM files + sample names  
./rnaedgeflow metagene \  
  --inputBam path/to/sample1.bam --sampleName Sample1 \  
  --inputBam2 path/to/sample2.bam --sampleName2 Sample2 \  
  --inputBam3 path/to/sample3.bam --sampleName3 Sample3 \  
  --outputDir results/metagene_out  

Options:  
--upstream <bp>    number of bases upstream of TSS (default: 1000)  
--downstream <bp>  number of bases downstream of TES (default: 1000)  
--bins <N>         number of bins per region (default: 100)  

Example:  
./rnaedgeflow metagene \  
  --inputDir example_data/test_001/process \  
  --upstream 500 --downstream 500 --bins 100 \  
  --outputDir example_data/test_001/result/metagene  

---

## Example Data

A small test dataset is included under `example_data/test_001`.  
You can test the pipeline with:  
./rnaedgeflow pipeline --profile example_data/example.profile  
./rnaedgeflow metagene --inputDir example_data/test_001/process --bins 100  

---

## Citation

If you use **RNAEdgeFlow** in your work, please cite this repository.  



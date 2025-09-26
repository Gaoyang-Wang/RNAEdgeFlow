# ================== Required parameters ==================
# Output directory for all results
OutputDir=./example_data

# Sample name (will be used as subdirectory under OutputDir)
SampleName=test_001

# Number of threads to use
Threads=10

# Length of UMI sequences
UMILength=8

# Directory containing raw FASTQ files
InputDir=./example_data/raw_fastq

# Prefix of raw FASTQ file names (before _R1/_R2)
RawFastqPre=BCnewAlldt_MaUP_YPS_S16

# Path to STAR genome index
ReferenceGenome=./genome/STARIndex2.7.11


# ================== 3' end (P3) ==================
# Barcode pattern for umi_tools (X = fixed base, N = UMI base)
P3_BC_Pattern=XXXXXXNNNNNNNNXXX

# UMI region (1-based coordinate, start:end)
P3_UmiRegion=7:14

# Fixed sequence for filtering (optional)
P3_FixedSeq1=
# File containing fixed sequence list (optional)
P3_FixedSeq1_File=./example_data/p3_barcode_position1_6.txt
P3_FixedSeq1Region=1:6
P3_FixedSeq1Mismatch=0

# Second fixed sequence (optional)
P3_FixedSeq2=TTT
P3_FixedSeq2_File=
P3_FixedSeq2Region=15:17
P3_FixedSeq2Mismatch=0


# ================== 5' end (P5) ==================
# Barcode pattern for umi_tools
P5_BC_Pattern=XXXXNNNNNNNNXXXXX
P5_UmiRegion=5:12

# Fixed sequence 1
P5_FixedSeq1=AGTC
P5_FixedSeq1Region=1:4
P5_FixedSeq1Mismatch=0

# Fixed sequence 2
P5_FixedSeq2=CATCAGGG
P5_FixedSeq2Region=13:20
P5_FixedSeq2Mismatch=0


# ================== Optional parameters ==================
# Minimum read length to keep after trimming
MinReadLen=5

# Extra STAR parameters (leave empty if not needed)
STAR_EXTRA=""

# Adapter trimming options for 3' reads (cutadapt)
CUTADAPT_3P_OPTS="--rc -m 5 -j 0 -n 4 \
-g ^CTGAGCTTTNNTTTNNTTTNN -g ^AGCAGCTTTNNTTTNNTTTNN -g ^GCTAGCTTTNNTTTNNTTTNN \
-g C{8} -G C{8} \
-a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

# Adapter trimming options for 5' reads
CUTADAPT_5P_OPTS="--rc -m 5 -j 0 -n 4 \
-g ^AGTCCATCAGGG -g C{8} \
-G C{8} \
-a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

# Adapter trimming options for internal reads
CUTADAPT_INTERNAL_OPTS="--rc -m 5 -j 0 -n 4 \
-g C{8} -a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-G C{8} -A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

# Number of bases to calculate base composition stats
BaseStatNumBases=85


# ================== Internal expression (optional: only if provided) ==================
# Gene annotation GTF (for StringTie)
GeneAnnotationGTF=./genome/sacCer3.ncbiRefSeq.gtf

# Path to prepDE.py (from StringTie package)
PrepDEPy=./software/prepDE.py

# Extra StringTie parameters (leave empty if not needed)
StringTieExtra=""

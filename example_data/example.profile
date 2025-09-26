# ================== 必填参数 ==================
OutputDir=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/project/EPIC/RNAEdgeFlow/example_data
SampleName=test_001
Threads=10
UMILength=8
InputDir=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/project/EPIC/RNAEdgeFlow/example_data/raw_fastq
RawFastqPre=BCnewAlldt_MaUP_YPS_S16
ReferenceGenome=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/project/EPIC/TSS_PAS/data/genome/STARIndex2.7.11

# ================== 3' 端 (P3) ==================
P3_BC_Pattern=XXXXXXNNNNNNNNXXX
P3_UmiRegion=7:14
P3_FixedSeq1=
P3_FixedSeq1_File=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/project/EPIC/RNAEdgeFlow/example_data/p3_barcode_position1_6.txt
P3_FixedSeq1Region=1:6
P3_FixedSeq1Mismatch=0
P3_FixedSeq2=TTT
P3_FixedSeq2_File=
P3_FixedSeq2Region=15:17
P3_FixedSeq2Mismatch=0

# ================== 5' 端 (P5) ==================
P5_BC_Pattern=XXXXNNNNNNNNXXXXX
P5_UmiRegion=5:12
P5_FixedSeq1=AGTC
P5_FixedSeq1Region=1:4
P5_FixedSeq1Mismatch=0
P5_FixedSeq2=CATCAGGG
P5_FixedSeq2Region=13:20
P5_FixedSeq2Mismatch=0

# ================== 选填参数 ==================
MinReadLen=5
STAR_EXTRA=""

CUTADAPT_3P_OPTS="--rc -m 5 -j 0 -n 4 \
-g ^CTGAGCTTTNNTTTNNTTTNN -g ^AGCAGCTTTNNTTTNNTTTNN -g ^GCTAGCTTTNNTTTNNTTTNN \
-g C{8} -G C{8} \
-a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

CUTADAPT_5P_OPTS="--rc -m 5 -j 0 -n 4 \
-g ^AGTCCATCAGGG -g C{8} \
-G C{8} \
-a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

CUTADAPT_INTERNAL_OPTS="--rc -m 5 -j 0 -n 4 \
-g C{8} -a G{8} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-G C{8} -A G{8} -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--quality-cutoff=20"

BaseStatNumBases=85

# Internal 表达量（填了才会启用）
GeneAnnotationGTF=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/project/EPIC/TSS_PAS/data/genome/sacCer3.ncbiRefSeq.gtf
PrepDEPy=/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/gywang/software/miniconda3/envs/readsMaster/bin/prepDE.py
StringTieExtra=""
#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Args: only --profile
# ---------------------------
PROFILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile=*) PROFILE="${1#*=}"; shift ;;
    --profile)    PROFILE="$2"; shift 2 ;;
    *) echo "[ERR] Unknown arg: $1" >&2; exit 1 ;;
  esac
done
[[ -n "${PROFILE}" && -f "${PROFILE}" ]] || { echo "[ERR] Profile missing: ${PROFILE}"; exit 1; }

# Load profile
# shellcheck source=/dev/null
source "${PROFILE}"

# ===========================
# Required keys in profile
# ===========================
req_vars=( OutputDir SampleName Threads UMILength InputDir RawFastqPre ReferenceGenome
           P3_BC_Pattern P3_UmiRegion P3_FixedSeq1Region P3_FixedSeq1Mismatch P3_FixedSeq2Region P3_FixedSeq2Mismatch
           P5_BC_Pattern P5_UmiRegion P5_FixedSeq1Region P5_FixedSeq1Mismatch P5_FixedSeq2Region P5_FixedSeq2Mismatch )
for v in "${req_vars[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "[ERR] Missing key in profile: $v"; exit 1; }
done

# Optional keys (defaults)
MinReadLen="${MinReadLen:-5}"
STAR_EXTRA="${STAR_EXTRA:-}"
CUTADAPT_3P_OPTS="${CUTADAPT_3P_OPTS:-}"
CUTADAPT_5P_OPTS="${CUTADAPT_5P_OPTS:-}"
CUTADAPT_INTERNAL_OPTS="${CUTADAPT_INTERNAL_OPTS:-}"
P3_FixedSeq1="${P3_FixedSeq1:-}"
P3_FixedSeq1_File="${P3_FixedSeq1_File:-}"
P3_FixedSeq2="${P3_FixedSeq2:-}"
P5_FixedSeq1="${P5_FixedSeq1:-}"
P5_FixedSeq2="${P5_FixedSeq2:-}"
# Internal expression (optional)
GeneAnnotationGTF="${GeneAnnotationGTF:-}"
PrepDEPy="${PrepDEPy:-}"
StringTieExtra="${StringTieExtra:-}"
# Base-content helper (optional; default 85)
BaseStatNumBases="${BaseStatNumBases:-85}"

THREADS="${Threads}"
sample="${SampleName}"

# ===========================
# Directory layout under OutputDir/SampleName
# ===========================
base_dir="${OutputDir}/${sample}"

process_dir="${base_dir}/process"
process_readsID_dir="${process_dir}/readsID"

result_dir="${base_dir}/result"
terminal_bed_dir="${result_dir}/terminal_bed"
terminal_bed_raw_dir="${terminal_bed_dir}/raw"
internal_expr_dir="${result_dir}/internal_expr"
stat_dir="${result_dir}/stat"

log_dir="${base_dir}/log"   # peer to result
mkdir -p "${process_dir}" "${process_readsID_dir}" \
         "${terminal_bed_raw_dir}" \
         "${internal_expr_dir}" "${stat_dir}" \
         "${log_dir}"

# --------- Logging: single PIPELOG only ---------
PIPELOG="${log_dir}/${sample}.pipeline.log"
: > "${PIPELOG}"

ts() { date +'%F %T'; }   # timestamp
CURRENT_STEP=""

fail_and_exit() {
  local msg="$1"
  echo "[FAIL] $(ts) [${CURRENT_STEP}] ${msg}" | tee -a "${PIPELOG}" >&2
  exit 1
}

on_error() {
  local ln="$1"; local cmd="$2"
  echo "[ERR ] $(ts) [${CURRENT_STEP}] line ${ln}: ${cmd}" >> "${PIPELOG}"
  echo "[FAIL] $(ts) [${CURRENT_STEP}] See ${PIPELOG}" >> "${PIPELOG}"
  exit 1
}
trap 'on_error $LINENO "$BASH_COMMAND"' ERR

step() {
  CURRENT_STEP="$1"
  echo "===== [$(ts)] STEP: ${CURRENT_STEP} =====" | tee -a "${PIPELOG}" >/dev/null
  echo "[START] $(ts) [${CURRENT_STEP}]" >> "${PIPELOG}"
}

ok() {
  local detail="${1:-}"
  echo "[OK  ] $(ts) [${CURRENT_STEP}] ${detail}" >> "${PIPELOG}"
}

check_file() {
  local f="$1"; local label="${2:-$1}"
  [[ -s "$f" ]] || fail_and_exit "Missing or empty file: ${label}"
}

# --------- Locate FASTQ with wildcards; require unique match ---------
pick_fastq() {
  local dir="$1" pre="$2" r12="$3"
  shopt -s nullglob
  local candidates=( \
    "${dir}/${pre}_${r12}"*.fastq.gz \
    "${dir}/${pre}_${r12}"*.fastq \
  )
  shopt -u nullglob
  local n="${#candidates[@]}"
  if (( n == 0 )); then
    return 1
  elif (( n > 1 )); then
    echo "[ERR ] Multiple ${r12} FASTQs matched for prefix '${pre}' in '${dir}':" >> "${PIPELOG}"
    printf ' - %s\n' "${candidates[@]}" >> "${PIPELOG}"
    echo "[HINT] Refine RawFastqPre or clean directory to ensure a single match." >> "${PIPELOG}"
    return 2
  else
    echo "${candidates[0]}"
    return 0
  fi
}

# ===========================
# 0) Locate input FASTQs and STAR index
# ===========================
step "Locate FASTQ and STAR index"
R1="$(pick_fastq "${InputDir}" "${RawFastqPre}" "R1")" || fail_and_exit "R1 not found or multiple matches"
R2="$(pick_fastq "${InputDir}" "${RawFastqPre}" "R2")" || fail_and_exit "R2 not found or multiple matches"
[[ -d "${ReferenceGenome}" ]] || fail_and_exit "STAR index not found: ${ReferenceGenome}"
echo "[INFO] Sample=${sample}" | tee -a "${PIPELOG}" >/dev/null
echo "[INFO] R1=${R1}" | tee -a "${PIPELOG}" >/dev/null
echo "[INFO] R2=${R2}" | tee -a "${PIPELOG}" >/dev/null
echo "[INFO] STAR index=${ReferenceGenome}" | tee -a "${PIPELOG}" >/dev/null
echo "[INFO] OutputDir=${OutputDir}" | tee -a "${PIPELOG}" >/dev/null
ok "FASTQ and index found"

# ===========================
# Helper: select R1 terminal reads using fixed motifs/positions
# ===========================
filter_terminal() {
  local which="$1" out_fastq="$2"
  local fixed1_file fixed1 fixed1_region fixed1_mm fixed2 fixed2_region fixed2_mm
  if [[ "${which}" == "3P" ]]; then
    fixed1_file="${P3_FixedSeq1_File}"
    fixed1="${P3_FixedSeq1}"
    fixed1_region="${P3_FixedSeq1Region}"
    fixed1_mm="${P3_FixedSeq1Mismatch}"
    fixed2="${P3_FixedSeq2}"
    fixed2_region="${P3_FixedSeq2Region}"
    fixed2_mm="${P3_FixedSeq2Mismatch}"
  else
    fixed1_file=""
    fixed1="${P5_FixedSeq1}"
    fixed1_region="${P5_FixedSeq1Region}"
    fixed1_mm="${P5_FixedSeq1Mismatch}"
    fixed2="${P5_FixedSeq2}"
    fixed2_region="${P5_FixedSeq2Region}"
    fixed2_mm="${P5_FixedSeq2Mismatch}"
  fi

  if [[ -n "${fixed1_file}" ]]; then
    seqkit grep -s -j "${THREADS}" -f "${fixed1_file}" -R "${fixed1_region}" -m "${fixed1_mm}" "${R1}"
  elif [[ -n "${fixed1}" ]]; then
    seqkit grep -s -j "${THREADS}" -p "${fixed1}" -R "${fixed1_region}" -m "${fixed1_mm}" "${R1}"
  else
    cat "${R1}"
  fi \
  | { if [[ -n "${fixed2}" ]]; then
        seqkit grep -s -j "${THREADS}" -p "${fixed2}" -R "${fixed2_region}" -m "${fixed2_mm}"
      else
        cat
      fi; } \
  | gzip -c > "${out_fastq}"
}

# ===========================
# 1) Extract terminal reads + split internal
# ===========================
step "Extract terminal reads (R1)"
f_3p_r1="${process_dir}/${sample}_3P_R1.fastq.gz"
f_5p_r1="${process_dir}/${sample}_5P_R1.fastq.gz"
filter_terminal "3P" "${f_3p_r1}" >> "${PIPELOG}" 2>&1
filter_terminal "5P" "${f_5p_r1}" >> "${PIPELOG}" 2>&1
check_file "${f_3p_r1}" "3P R1"; check_file "${f_5p_r1}" "5P R1"
ok "3P/5P R1 extracted"

step "Collect read IDs and split internal (R1/R2)"
id_3p_r1="${process_readsID_dir}/${sample}_3P_R1_readsID.txt"
id_5p_r1="${process_readsID_dir}/${sample}_5P_R1_readsID.txt"
seqkit seq -n "${f_3p_r1}" > "${id_3p_r1}" 2>> "${PIPELOG}"
seqkit seq -n "${f_5p_r1}" > "${id_5p_r1}" 2>> "${PIPELOG}"
check_file "${id_3p_r1}" "3P R1 IDs"; check_file "${id_5p_r1}" "5P R1 IDs"
all_term_r1="${process_readsID_dir}/${sample}_R1_terminal_readsID.txt"
cat "${id_3p_r1}" "${id_5p_r1}" > "${all_term_r1}"
check_file "${all_term_r1}" "R1 terminal IDs"

f_int_r1="${process_dir}/${sample}_internal_R1.fastq.gz"
seqkit grep -n -j "${THREADS}" -v -f "${all_term_r1}" "${R1}" -o "${f_int_r1}" >> "${PIPELOG}" 2>&1
check_file "${f_int_r1}" "internal R1"

id_3p_r2="${process_readsID_dir}/${sample}_3P_R2_readsID.txt"
id_5p_r2="${process_readsID_dir}/${sample}_5P_R2_readsID.txt"
all_term_r2="${process_readsID_dir}/${sample}_R2_terminal_readsID.txt"
sed 's/1:N/2:N/' "${id_3p_r1}" > "${id_3p_r2}" 2>> "${PIPELOG}"
sed 's/1:N/2:N/' "${id_5p_r1}" > "${id_5p_r2}" 2>> "${PIPELOG}"
sed 's/1:N/2:N/' "${all_term_r1}" > "${all_term_r2}" 2>> "${PIPELOG}"
check_file "${id_3p_r2}" "3P R2 IDs"; check_file "${id_5p_r2}" "5P R2 IDs"; check_file "${all_term_r2}" "R2 terminal IDs"

f_3p_r2="${process_dir}/${sample}_3P_R2.fastq.gz"
f_5p_r2="${process_dir}/${sample}_5P_R2.fastq.gz"
f_int_r2="${process_dir}/${sample}_internal_R2.fastq.gz"
seqkit grep -n -j "${THREADS}" -f "${id_3p_r2}" "${R2}" -o "${f_3p_r2}" >> "${PIPELOG}" 2>&1
seqkit grep -n -j "${THREADS}" -f "${id_5p_r2}" "${R2}" -o "${f_5p_r2}" >> "${PIPELOG}" 2>&1
seqkit grep -n -j "${THREADS}" -v -f "${all_term_r2}" "${R2}" -o "${f_int_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_3p_r2}" "3P R2"; check_file "${f_5p_r2}" "5P R2"; check_file "${f_int_r2}" "internal R2"
ok "IDs exported and internal/terminal split (R1/R2)"

rm -f "${id_3p_r1}" "${id_5p_r1}" "${id_3p_r2}" "${id_5p_r2}" "${all_term_r1}" "${all_term_r2}"

# ===========================
# 1.5) Base content & read-length stats for RAW R1 (3P/5P/internal)
# ===========================
step "Base content/length stats (raw R1: 3P/5P/internal)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_STAT_PY="${SCRIPT_DIR}/base_content_and_read_length_stat.sh"
if [[ ! -f "${BASE_STAT_PY}" ]]; then
  echo "[WARN] base_content_and_read_length_stat.sh not found at ${BASE_STAT_PY}; skip this step." | tee -a "${PIPELOG}"
  ok "Skipped (script not found)"
else
  run_base_stat() {
    local in_fq="$1" tag="$2"
    if [[ -s "${in_fq}" ]]; then
      echo "[RUN ] base_stat ${tag}: ${in_fq}" >> "${PIPELOG}"
      python "${BASE_STAT_PY}" \
        -i "${in_fq}" \
        -n "${BaseStatNumBases}" \
        -o "${stat_dir}" \
        >> "${PIPELOG}" 2>&1
      echo "[OK  ] base_stat ${tag}" >> "${PIPELOG}"
    else
      echo "[WARN] Skip base_stat ${tag}: missing FASTQ ${in_fq}" >> "${PIPELOG}"
    fi
  }
  run_base_stat "${f_3p_r1}" "3P_R1_raw"
  run_base_stat "${f_5p_r1}" "5P_R1_raw"
  run_base_stat "${f_int_r1}" "internal_R1_raw"
  ok "Base stats generated (see ${stat_dir})"
fi

# ===========================
# 2) UMI extraction
# ===========================
step "UMI extraction"
f_3p_umi_r1="${process_dir}/${sample}_3P_umi_R1.fastq.gz"
f_3p_umi_r2="${process_dir}/${sample}_3P_umi_R2.fastq.gz"
umi_tools extract --bc-pattern="${P3_BC_Pattern}" -I "${f_3p_r1}" --read2-in "${f_3p_r2}" --stdout "${f_3p_umi_r1}" --read2-out "${f_3p_umi_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_3p_umi_r1}" "3P umi R1"; check_file "${f_3p_umi_r2}" "3P umi R2"

f_5p_umi_r1="${process_dir}/${sample}_5P_umi_R1.fastq.gz"
f_5p_umi_r2="${process_dir}/${sample}_5P_umi_R2.fastq.gz"
umi_tools extract --bc-pattern="${P5_BC_Pattern}" -I "${f_5p_r1}" --read2-in "${f_5p_r2}" --stdout "${f_5p_umi_r1}" --read2-out "${f_5p_umi_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_5p_umi_r1}" "5P umi R1"; check_file "${f_5p_umi_r2}" "5P umi R2"
ok "UMI extracted"

# ===========================
# 3) Adapter trimming + fastp
# ===========================
step "Adapter trimming + fastp"
f_3p_cut_r1="${process_dir}/${sample}_3P_umi_cutted_R1.fastq.gz"
f_3p_cut_r2="${process_dir}/${sample}_3P_umi_cutted_R2.fastq.gz"
cutadapt ${CUTADAPT_3P_OPTS} -o "${f_3p_cut_r1}" -p "${f_3p_cut_r2}" "${f_3p_umi_r1}" "${f_3p_umi_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_3p_cut_r1}" "3P cut R1"; check_file "${f_3p_cut_r2}" "3P cut R2"
f_3p_clean_r1="${process_dir}/${sample}_3P_clean_R1.fastq.gz"
f_3p_clean_r2="${process_dir}/${sample}_3P_clean_R2.fastq.gz"
fastp --thread "${THREADS}" --length_required "${MinReadLen}" -i "${f_3p_cut_r1}" -I "${f_3p_cut_r2}" -o "${f_3p_clean_r1}" -O "${f_3p_clean_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_3p_clean_r1}" "3P clean R1"; check_file "${f_3p_clean_r2}" "3P clean R2"

f_5p_cut_r1="${process_dir}/${sample}_5P_umi_cutted_R1.fastq.gz"
f_5p_cut_r2="${process_dir}/${sample}_5P_umi_cutted_R2.fastq.gz"
cutadapt ${CUTADAPT_5P_OPTS} -o "${f_5p_cut_r1}" -p "${f_5p_cut_r2}" "${f_5p_umi_r1}" "${f_5p_umi_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_5p_cut_r1}" "5P cut R1"; check_file "${f_5p_cut_r2}" "5P cut R2"
f_5p_clean_r1="${process_dir}/${sample}_5P_clean_R1.fastq.gz"
f_5p_clean_r2="${process_dir}/${sample}_5P_clean_R2.fastq.gz"
fastp --thread "${THREADS}" --length_required "${MinReadLen}" -i "${f_5p_cut_r1}" -I "${f_5p_cut_r2}" -o "${f_5p_clean_r1}" -O "${f_5p_clean_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_5p_clean_r1}" "5P clean R1"; check_file "${f_5p_clean_r2}" "5P clean R2"

f_int_cut_r1="${process_dir}/${sample}_internal_cutted_R1.fastq.gz"
f_int_cut_r2="${process_dir}/${sample}_internal_cutted_R2.fastq.gz"
cutadapt ${CUTADAPT_INTERNAL_OPTS} -o "${f_int_cut_r1}" -p "${f_int_cut_r2}" "${f_int_r1}" "${f_int_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_int_cut_r1}" "internal cut R1"; check_file "${f_int_cut_r2}" "internal cut R2"
f_int_clean_r1="${process_dir}/${sample}_internal_clean_R1.fastq.gz"
f_int_clean_r2="${process_dir}/${sample}_internal_clean_R2.fastq.gz"
fastp --thread "${THREADS}" --length_required "${MinReadLen}" -i "${f_int_cut_r1}" -I "${f_int_cut_r2}" -o "${f_int_clean_r1}" -O "${f_int_clean_r2}" >> "${PIPELOG}" 2>&1
check_file "${f_int_clean_r1}" "internal clean R1"; check_file "${f_int_clean_r2}" "internal clean R2"
ok "cutadapt + fastp done"

# Cleanup intermediates
rm -f "${f_3p_umi_r1}" "${f_3p_umi_r2}" "${f_5p_umi_r1}" "${f_5p_umi_r2}" \
      "${f_3p_cut_r1}" "${f_3p_cut_r2}" "${f_5p_cut_r1}" "${f_5p_cut_r2}" \
      "${f_int_cut_r1}" "${f_int_cut_r2}"

# ===========================
# 4) STAR mapping + UMI dedup for 3P/5P
# ===========================
step "STAR alignment"
run_map() {
  local tag="$1"
  local r1="$2"
  local r2="$3"

  local prefix="${process_dir}/${sample}_${tag}_star_"
  local outbam="${prefix}Aligned.sortedByCoord.out.bam"

  local READCMD="cat"
  [[ "${r1}" == *.gz ]] && READCMD="zcat"

  STAR --runThreadN "${THREADS}" \
       --runMode alignReads \
       --readFilesCommand "${READCMD}" \
       --genomeDir "${ReferenceGenome}" \
       --genomeLoad LoadAndKeep \
       --readFilesIn "${r1}" "${r2}" \
       --limitBAMsortRAM 100000000000 \
       --alignIntronMax 2500 \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMultimapNmax 1 \
       --outFileNamePrefix "${prefix}" \
       --outStd Log \
       --outBAMsortingThreadN 6 \
       ${STAR_EXTRA} \
       >> "${PIPELOG}" 2>&1

  check_file "${outbam}" "${tag} STAR BAM"
  samtools index --threads "${THREADS}" "${outbam}" >> "${PIPELOG}" 2>&1
}
run_map "3P" "${f_3p_clean_r1}" "${f_3p_clean_r2}"
run_map "5P" "${f_5p_clean_r1}" "${f_5p_clean_r2}"
run_map "internal" "${f_int_clean_r1}" "${f_int_clean_r2}"
ok "STAR BAMs generated"

step "UMI dedup (3P/5P) and finalize BAMs"
for tag in 3P 5P; do
  inbam="${process_dir}/${sample}_${tag}_star_Aligned.sortedByCoord.out.bam"
  dedup="${process_dir}/${sample}_${tag}_star_dedup.bam"
  umi_tools dedup -I "${inbam}" -S "${dedup}" --method=unique --paired >> "${PIPELOG}" 2>&1
  check_file "${dedup}" "${tag} dedup BAM"
  rm -f "${inbam}" "${inbam}.bai"
  mv "${dedup}" "${process_dir}/${sample}_${tag}_final.bam"
  samtools index --threads "${THREADS}" "${process_dir}/${sample}_${tag}_final.bam" >> "${PIPELOG}" 2>&1
  check_file "${process_dir}/${sample}_${tag}_final.bam.bai" "${tag} final BAM index"
done
mv "${process_dir}/${sample}_internal_star_Aligned.sortedByCoord.out.bam" \
   "${process_dir}/${sample}_internal_final.bam"
samtools index --threads "${THREADS}" "${process_dir}/${sample}_internal_final.bam" >> "${PIPELOG}" 2>&1
check_file "${process_dir}/${sample}_internal_final.bam.bai" "internal final BAM index"
ok "Final BAMs ready"

for tag in 3P 5P internal; do
  samtools flagstat --threads "${THREADS}" \
    "${process_dir}/${sample}_${tag}_final.bam" \
    > "${stat_dir}/${sample}_${tag}.flagstat" 2>> "${PIPELOG}"
  check_file "${stat_dir}/${sample}_${tag}.flagstat" "${tag} flagstat"
done

# ===========================
# X) Single-sample read counts & fractions
# ===========================
step "Sample read statistics (counts & fractions only)"

percent() {
  if [[ "${2:-0}" -eq 0 ]]; then
    echo "0.00"
  else
    echo "scale=3; ${1}/${2}" | bc | awk '{printf "%.2f\n", $1}'
  fi
}

count_reads_only() {
  local fq="$1"
  if [[ -f "$fq" ]]; then
    seqkit stats -j "${THREADS}" -a "$fq" | awk 'NR==2 {gsub(",","",$4); print $4}'
  else
    echo "0"
  fi
}

count_mapped_reads() {
  local bam="$1"
  if [[ -f "$bam" ]]; then
    samtools flagstat --threads "${THREADS}" "$bam" | awk 'NR==7 {print $1}'
  else
    echo "0"
  fi
}

# raw (R1)
raw_all_n="$(count_reads_only "${R1}")"
raw_3p_n="$(count_reads_only "${f_3p_r1}")"
raw_5p_n="$(count_reads_only "${f_5p_r1}")"
raw_int_n="$(count_reads_only "${f_int_r1}")"

# clean (R1)
clean_3p_n="$(count_reads_only "${f_3p_clean_r1}")"
clean_5p_n="$(count_reads_only "${f_5p_clean_r1}")"
clean_int_n="$(count_reads_only "${f_int_clean_r1}")"
clean_all_n=$(( clean_3p_n + clean_5p_n + clean_int_n ))

# mapped (BAM)
mapped_3p_n="$(count_mapped_reads "${process_dir}/${sample}_3P_final.bam")"
mapped_5p_n="$(count_mapped_reads "${process_dir}/${sample}_5P_final.bam")"
mapped_int_n="$(count_mapped_reads "${process_dir}/${sample}_internal_final.bam")"
mapped_all_n=$(( mapped_3p_n + mapped_5p_n + mapped_int_n ))

# fractions
raw_3p_pct="$(percent "${raw_3p_n}"   "${raw_all_n}")"
raw_5p_pct="$(percent "${raw_5p_n}"   "${raw_all_n}")"
raw_int_pct="$(percent "${raw_int_n}" "${raw_all_n}")"

clean_3p_pct="$(percent "${clean_3p_n}"   "${clean_all_n}")"
clean_5p_pct="$(percent "${clean_5p_n}"   "${clean_all_n}")"
clean_int_pct="$(percent "${clean_int_n}" "${clean_all_n}")"

mapped_3p_pct="$(percent "${mapped_3p_n}"   "${mapped_all_n}")"
mapped_5p_pct="$(percent "${mapped_5p_n}"   "${mapped_all_n}")"
mapped_int_pct="$(percent "${mapped_int_n}" "${mapped_all_n}")"

reads_stat_tsv="${stat_dir}/${sample}_reads_count.tsv"
{
  echo -e "sample\tpart\tstage\treads\tfraction"
  echo -e "${sample}\tall\traw\t${raw_all_n}\t1.00"
  echo -e "${sample}\t3P\traw\t${raw_3p_n}\t${raw_3p_pct}"
  echo -e "${sample}\t5P\traw\t${raw_5p_n}\t${raw_5p_pct}"
  echo -e "${sample}\tinternal\traw\t${raw_int_n}\t${raw_int_pct}"

  echo -e "${sample}\tall\tclean\t${clean_all_n}\t1.00"
  echo -e "${sample}\t3P\tclean\t${clean_3p_n}\t${clean_3p_pct}"
  echo -e "${sample}\t5P\tclean\t${clean_5p_n}\t${clean_5p_pct}"
  echo -e "${sample}\tinternal\tclean\t${clean_int_n}\t${clean_int_pct}"

  echo -e "${sample}\tall\tmapped\t${mapped_all_n}\t1.00"
  echo -e "${sample}\t3P\tmapped\t${mapped_3p_n}\t${mapped_3p_pct}"
  echo -e "${sample}\t5P\tmapped\t${mapped_5p_n}\t${mapped_5p_pct}"
  echo -e "${sample}\tinternal\tmapped\t${mapped_int_n}\t${mapped_int_pct}"
} > "${reads_stat_tsv}" 2>> "${PIPELOG}"

check_file "${reads_stat_tsv}" "reads_count.tsv"
ok "Reads stats saved -> ${reads_stat_tsv}"

# ===========================
# 5) Coverage bedgraphs (raw) -> result/terminal_bed
# ===========================
step "Coverage (raw only)"
cov_do() {
  local tag="$1" bam="${process_dir}/${sample}_${tag}_final.bam"
  bedtools genomecov -ibam "${bam}" -5 -bg \
    > "${terminal_bed_raw_dir}/${sample}_${tag}.bedgraph" 2>> "${PIPELOG}"
  bedtools genomecov -ibam "${bam}" -5 -strand + -bg \
    > "${terminal_bed_raw_dir}/${sample}_${tag}_plus.bedgraph" 2>> "${PIPELOG}"
  bedtools genomecov -ibam "${bam}" -5 -strand - -bg \
    > "${terminal_bed_raw_dir}/${sample}_${tag}_minus.bedgraph" 2>> "${PIPELOG}"

  check_file "${terminal_bed_raw_dir}/${sample}_${tag}.bedgraph"       "${tag} raw cov"
  check_file "${terminal_bed_raw_dir}/${sample}_${tag}_plus.bedgraph"  "${tag} raw plus"
  check_file "${terminal_bed_raw_dir}/${sample}_${tag}_minus.bedgraph" "${tag} raw minus"
}
cov_do 5P
cov_do 3P
ok "Coverage generated (raw)"

# ===========================
# 6) Internal expression (StringTie + prepDE)
# ===========================
step "Internal expression (StringTie + prepDE)"
if [[ -z "${GeneAnnotationGTF}" || -z "${PrepDEPy}" ]]; then
  echo "[SKIP] GeneAnnotationGTF or PrepDEPy not set in profile; skip internal expression." | tee -a "${PIPELOG}"
  ok "Skipped (no GTF or prepDE.py)"
else
  bam_int="${process_dir}/${sample}_internal_final.bam"
  check_file "${bam_int}" "internal final BAM"

  expr_tsv="${internal_expr_dir}/${sample}_internal_geneExp.tsv"
  sample_gtf="${internal_expr_dir}/${sample}_internal_sample.gtf"
  listfile="${internal_expr_dir}/stringtie_samples.txt"

  stringtie "${bam_int}" -B -p "${THREADS}" -e \
    -G "${GeneAnnotationGTF}" -A "${expr_tsv}" \
    -o "${sample_gtf}" ${StringTieExtra} \
    >> "${PIPELOG}" 2>&1

  check_file "${expr_tsv}"  "stringtie -A tsv"
  check_file "${sample_gtf}" "stringtie sample gtf"

  touch "${listfile}"
  grep -v -F -e "${sample}" "${listfile}" > "${listfile}.tmp" || true
  mv "${listfile}.tmp" "${listfile}"
  printf "%s\t%s\n" "${sample}" "${sample_gtf}" >> "${listfile}"

  gene_mat="${internal_expr_dir}/gene_count_matrix.csv"
  tx_mat="${internal_expr_dir}/transcript_count_matrix.csv"
  python "${PrepDEPy}" -i "${listfile}" -g "${gene_mat}" -t "${tx_mat}" >> "${PIPELOG}" 2>&1

  check_file "${gene_mat}" "gene_count_matrix.csv"
  check_file "${tx_mat}"   "transcript_count_matrix.csv"
  ok "Internal expression matrices ready"
fi

# ===========================
# 7) QC (fastqc/multiqc) -> result/stat
# ===========================
step "QC (fastqc/multiqc)"
fastqc --threads "${THREADS}" -o "${stat_dir}" "${R1}" "${R2}" >> "${PIPELOG}" 2>&1
multiqc -f -o "${stat_dir}" "${stat_dir}" --filename "${sample}_multiqc" >> "${PIPELOG}" 2>&1
if compgen -G "${stat_dir}/*fastqc.html" > /dev/null; then :; else
  fail_and_exit "No fastqc html found in ${stat_dir}"
fi
check_file "${stat_dir}/${sample}_multiqc.html" "multiqc html"
ok "QC reports ready"

echo "[DONE] ${sample}" | tee -a "${PIPELOG}" >/dev/null
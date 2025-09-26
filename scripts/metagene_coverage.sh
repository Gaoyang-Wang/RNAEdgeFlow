#!/usr/bin/env bash
set -euo pipefail

ts(){ date +'%F %T'; }
log(){ echo "[$(ts)] $*"; }

# ---------- args ----------
PROFILE=""
GTF=""
INPUTDIR=""
OUTDIR=""            
BINS=100
UP=1000
DOWN=1000
declare -a BAMS=()
declare -a NAMES=()

while (( "$#" )); do
  case "$1" in
    --profile) PROFILE="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --inputdir) INPUTDIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;         
    --inputBam1|--inputBam2|--inputBam3) BAMS+=("$2"); shift 2;;
    --sample1|--sample2|--sample3) NAMES+=("$2"); shift 2;;
    --bins) BINS="$2"; shift 2;;
    --up)   UP="$2";   shift 2;;
    --down) DOWN="$2"; shift 2;;
    *)
      echo "[ERR] Unknown opt: $1" >&2
      exit 1;;
  esac
done

# ---------- resolve from profile (if given) ----------
if [[ -n "${PROFILE}" ]]; then
  # shellcheck disable=SC1090
  source "${PROFILE}"
  : "${OutputDir:?}"; : "${SampleName:?}"; : "${GeneAnnotationGTF:?}"

  GTF="${GTF:-$GeneAnnotationGTF}"
  out_root="${OutputDir%/}/${SampleName}"
  proc_dir="${out_root}/process"
  res_dir="${out_root}/result/metagene"
  log_dir="${out_root}/log"
  mkdir -p "${res_dir}/bigwig" "${log_dir}"

  auto_bams=(
    "${proc_dir}/${SampleName}_3P_final.bam"
    "${proc_dir}/${SampleName}_5P_final.bam"
    "${proc_dir}/${SampleName}_internal_final.bam"
  )
  auto_names=("3P" "5P" "internal")
  BAMS=(); NAMES=()
  for i in "${!auto_bams[@]}"; do
    if [[ -s "${auto_bams[$i]}" ]]; then
      BAMS+=("${auto_bams[$i]}")
      NAMES+=("${auto_names[$i]}")
    fi
  done

else
  if [[ -n "${OUTDIR}" ]]; then
    res_dir="${OUTDIR%/}"
  else
    res_dir="$(pwd)/metagene_out"
  fi
  log_dir="${res_dir}"
  mkdir -p "${res_dir}/bigwig" "${log_dir}"
fi

# ---------- alternative: inputdir ----------
if [[ -n "${INPUTDIR}" ]]; then
  mapfile -t found < <(find "${INPUTDIR}" -maxdepth 1 -type f -name "*.bam" | sort)
  if [[ "${#found[@]}" -eq 0 ]]; then
    echo "[ERR] No BAMs under ${INPUTDIR}" >&2; exit 1
  fi
  BAMS=("${found[@]}")
  NAMES=()
  for b in "${BAMS[@]}"; do
    bn="$(basename "${b}")"
    NAMES+=("${bn%.bam}")
  done
fi

# ---------- sanity ----------
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "[ERR] No BAMs provided/found. Use --profile or --inputdir or --inputBam*/--sample*." >&2
  exit 1
fi
if [[ -z "${GTF}" || ! -s "${GTF}" ]]; then
  echo "[ERR] Need --gtf <file> or set GeneAnnotationGTF in profile." >&2
  exit 1
fi

# ---------- R deps check/install ----------
Rscript - <<'RS' || true
pkgs <- c("optparse","GenomicRanges","GenomicFeatures","rtracklayer",
          "genomation","ggplot2","cowplot","txdbmaker","BiocManager")
need <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(need)>0) {
  if (!"BiocManager" %in% rownames(installed.packages()))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  bioc <- c("GenomicRanges","GenomicFeatures","rtracklayer","genomation","txdbmaker")
  tod_bioc <- intersect(need, bioc)
  tod_cran <- setdiff(need, union("BiocManager", bioc))
  if (length(tod_bioc)>0) BiocManager::install(tod_bioc, ask=FALSE, update=FALSE)
  if (length(tod_cran)>0) install.packages(tod_cran, repos="https://cloud.r-project.org")
}
RS

# ---------- chrom.sizes ----------
cs="${res_dir}/chrom.sizes"
if [[ ! -s "${cs}" ]]; then
  samtools view -H "${BAMS[0]}" \
  | awk '$1=="@SQ"{for(i=1;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4); if($i~/^LN:/)ln=substr($i,4)}; print sn"\t"ln}' > "${cs}"
fi

# ---------- per BAM: bedgraph(+/-) -> bigwig ----------
tbl="${res_dir}/samples.tsv"
echo -e "sample\tbwPlus\tbwMinus" > "${tbl}"

for i in "${!BAMS[@]}"; do
  bam="${BAMS[$i]}"; name="${NAMES[$i]}"
  base="${res_dir}/bigwig/${name}"
  bgp="${base}.plus.bedgraph"; bgm="${base}.minus.bedgraph"
  bwp="${base}.plus.bw";      bwm="${base}.minus.bw"

  log "[STEP] ${name}: BAM -> bigWig"
  bedtools genomecov -ibam "${bam}" -5 -strand + -bg > "${bgp}"
  bedtools genomecov -ibam "${bam}" -5 -strand - -bg > "${bgm}"
  bedGraphToBigWig "${bgp}" "${cs}" "${bwp}"
  bedGraphToBigWig "${bgm}" "${cs}" "${bwm}"
  echo -e "${name}\t${bwp}\t${bwm}" >> "${tbl}"
done

# ---------- call R (multi-sample TSV) ----------
png="${res_dir}/metagene_bw.png"
Rscript "$(dirname "$0")/metagene_plot_bw.R" \
  --gtf "${GTF}" --table "${tbl}" \
  --up "${UP}" --down "${DOWN}" --bins "${BINS}" \
  --out "${png}"

log "[DONE] Plot -> ${png}"
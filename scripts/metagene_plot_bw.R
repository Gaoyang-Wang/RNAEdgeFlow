#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(genomation)
  library(ggplot2); theme_set(theme_minimal())
})

optlist <- list(
  make_option("--gtf",       type="character"),
  make_option("--up",        type="integer",   default=1000),
  make_option("--down",      type="integer",   default=1000),
  make_option("--bins",      type="integer",   default=100),
  make_option("--out",       type="character", help="output PNG"),
  # legacy single-3tag mode 
  make_option("--bw3pPlus",  type="character"),
  make_option("--bw3pMinus", type="character"),
  make_option("--bw5pPlus",  type="character"),
  make_option("--bw5pMinus", type="character"),
  make_option("--bwintPlus", type="character"),
  make_option("--bwintMinus",type="character"),
  # new: multi-sample via TSV
  make_option("--table",     type="character", help="TSV with columns: sample,bwPlus,bwMinus")
)
opt <- parse_args(OptionParser(option_list = optlist))

if (!nzchar(opt$gtf)) stop("Need --gtf")
if (!nzchar(opt$out)) stop("Need --out")

# TxDb (Bioc 3.19+ uses txdbmaker)
txdb  <- txdbmaker::makeTxDbFromGFF(opt$gtf)
genes <- GenomicFeatures::genes(txdb)

nb <- as.integer(opt$bins)

load_bw_pair <- function(plus_bw, minus_bw) {
  gp <- rtracklayer::import(plus_bw,  format="BigWig")
  gm <- rtracklayer::import(minus_bw, format="BigWig")
  gp <- GenomicRanges::resize(gp, width=1, fix="start")
  gm <- GenomicRanges::resize(gm, width=1, fix="start")
  strand(gp) <- "+"
  strand(gm) <- "-"
  gr <- c(gp, gm)
  if (!"score" %in% names(mcols(gr))) mcols(gr)$score <- mcols(gr)[,1]
  gr
}

score_three_parts <- function(gr, genes, up_bp, down_bp, bins_each) {
  up_windows   <- GenomicRanges::promoters(genes, upstream=up_bp, downstream=0)
  body_windows <- genes
  dn_windows   <- GenomicRanges::flank(genes, width=down_bp, start=FALSE)

  up_mat <- genomation::ScoreMatrixBin(target=gr, windows=up_windows,   bin.num=bins_each, weight.col="score", strand.aware=TRUE)
  bd_mat <- genomation::ScoreMatrixBin(target=gr, windows=body_windows, bin.num=bins_each, weight.col="score", strand.aware=TRUE)
  dn_mat <- genomation::ScoreMatrixBin(target=gr, windows=dn_windows,   bin.num=bins_each, weight.col="score", strand.aware=TRUE)

  c(colMeans(up_mat), colMeans(bd_mat), colMeans(dn_mat))
}

scale_q95 <- function(v) {
  q <- as.numeric(stats::quantile(v, 0.95, na.rm = TRUE))
  if (!is.finite(q) || q <= 0) q <- max(v, na.rm = TRUE)
  if (!is.finite(q) || q <= 0) q <- 1
  v / q
}

# --------- build data.frame (multi or legacy) ----------
df <- NULL
x  <- seq_len(3*nb)

if (nzchar(opt$table)) {
  tab <- read.delim(opt$table, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  need <- c("sample","bwPlus","bwMinus")
  if (!all(need %in% names(tab))) stop("--table must have columns: sample,bwPlus,bwMinus")

  dd <- lapply(seq_len(nrow(tab)), function(i){
    smp <- tab$sample[i]
    gr  <- load_bw_pair(tab$bwPlus[i], tab$bwMinus[i])
    val <- score_three_parts(gr, genes, opt$up, opt$down, nb)
    val <- scale_q95(val)
    data.frame(x=x, value=val, tag=smp, stringsAsFactors=FALSE)
  })
  df <- do.call(rbind, dd)

} else {
  # legacy 3 curves
  stopifnot(nzchar(opt$bw3pPlus),  nzchar(opt$bw3pMinus),
            nzchar(opt$bw5pPlus),  nzchar(opt$bw5pMinus),
            nzchar(opt$bwintPlus), nzchar(opt$bwintMinus))
  gr_3p  <- load_bw_pair(opt$bw3pPlus,  opt$bw3pMinus)
  gr_5p  <- load_bw_pair(opt$bw5pPlus,  opt$bw5pMinus)
  gr_int <- load_bw_pair(opt$bwintPlus, opt$bwintMinus)
  v_3p   <- scale_q95(score_three_parts(gr_3p,  genes, opt$up, opt$down, nb))
  v_5p   <- scale_q95(score_three_parts(gr_5p,  genes, opt$up, opt$down, nb))
  v_int  <- scale_q95(score_three_parts(gr_int, genes, opt$up, opt$down, nb))
  df <- rbind(
    data.frame(x=x, value=v_5p,  tag="5P"),
    data.frame(x=x, value=v_int, tag="internal"),
    data.frame(x=x, value=v_3p,  tag="3P")
  )
  df$tag <- factor(df$tag, levels=c("5P","internal","3P"))
}

# --------- plot single panel ---------
p <- ggplot(df, aes(x=x, y=value, color=tag)) +
  geom_line(size=0.9) +
  geom_vline(xintercept = c(nb, 2*nb), linetype="dashed") +
  scale_x_continuous(
    breaks = c(1, nb, 2*nb, 3*nb),
    labels = c("Up start", "TSS", "Body|Down", "Down end")
  ) +
  labs(x = "Relative position (Up | Body | Down; equal bins)",
       y = "Normalized mean coverage (per-curve q95)",
       color = "Sample") +
  theme(legend.position = "right")

ggsave(filename = opt$out, plot = p, width = 10, height = 5, dpi = 180)
message("Saved plot -> ", opt$out)
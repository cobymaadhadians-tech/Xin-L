#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages(library(data.table))
records <- readRDS(file.path(root, "results/p5/neurological_phbc/formal/neurological_phbc_records.rds"))
n_block <- 200L
left_key <- "AD_GCST90704647_NOPROXY"
right_key <- "AD_GCST90704648_NOBIOBANK"
left <- records[[left_key]]
right <- records[[right_key]]
stopifnot(left$n_block == n_block, right$n_block == n_block,
          identical(left$auxiliaries, right$auxiliaries),
          length(left$block_uncorrected) == n_block,
          length(right$block_uncorrected) == n_block)

delta <- left$phbc - right$phbc
delta_jk <- left$corrected_weight^2 * left$block_uncorrected -
  right$corrected_weight^2 * right$block_uncorrected
pseudo <- n_block * delta - (n_block - 1L) * delta_jk
se <- sd(pseudo) / sqrt(n_block)
z <- delta / se
row <- data.table(
  comparison = "AD_NOPROXY_minus_NOBIOBANK",
  analysis_tier = "neurological_sensitivity_triangle_completion",
  left_target = left_key,
  right_target = right_key,
  left_phbc = left$phbc,
  right_phbc = right$phbc,
  delta_pp = 100 * delta,
  paired_jackknife_se_pp = 100 * se,
  ci95_low_pp = 100 * (delta - 1.96 * se),
  ci95_high_pp = 100 * (delta + 1.96 * se),
  z = z,
  p_two_sided = 2 * pnorm(-abs(z)),
  n_block = n_block,
  left_corrected_weight = left$corrected_weight,
  right_corrected_weight = right$corrected_weight,
  method = "200-block paired delete-one-block jackknife with target-specific fixed full-data correction weights"
)
out_dir <- file.path(root, "results/p6/supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(row, file.path(out_dir, "ad_noproxy_minus_nobiobank_paired_phbc.tsv"), sep = "\t")
cat("P6_AD_THIRD_PAIR_COMPLETE\n")

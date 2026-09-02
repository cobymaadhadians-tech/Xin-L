#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages(library(data.table))
records <- readRDS(file.path(root, "results/p25/final/p25_phbc_records.rds"))$full
n_block <- 200L
left_key <- "SCZ_FINNGEN_R13"
right_key <- "SCZ_PGC2022"
left <- records[[left_key]]
right <- records[[right_key]]
stopifnot(identical(left$auxiliaries, right$auxiliaries),
          length(left$block_uncorrected) == n_block,
          length(right$block_uncorrected) == n_block)

# Correction-off estimates use the uncorrected full-data quadratic forms and
# their paired leave-one-block-out analogues.
left_point <- left$uncorrected_quadratic_form
right_point <- right$uncorrected_quadratic_form
delta <- left_point - right_point
delta_jk <- left$block_uncorrected - right$block_uncorrected
pseudo <- n_block * delta - (n_block - 1L) * delta_jk
se <- sd(pseudo) / sqrt(n_block)
corrected_delta <- left$phbc - right$phbc
row <- data.table(
  comparison = "SCZ_FINNGEN_minus_PGC_correction_off",
  analysis_tier = "psychiatric_sensitivity",
  auxiliaries = paste(left$auxiliaries, collapse = ","),
  n_aux = length(left$auxiliaries),
  left_uncorrected_phbc = left_point,
  right_uncorrected_phbc = right_point,
  delta_uncorrected_pp = 100 * delta,
  paired_jackknife_se_pp = 100 * se,
  ci95_low_pp = 100 * (delta - 1.96 * se),
  ci95_high_pp = 100 * (delta + 1.96 * se),
  z = delta / se,
  p_two_sided = 2 * pnorm(-abs(delta / se)),
  corrected_delta_pp = 100 * corrected_delta,
  left_corrected_weight = left$corrected_weight,
  right_corrected_weight = right$corrected_weight,
  n_block = n_block,
  method = "uncorrected full-data quadratic form with 200-block paired delete-one-block jackknife"
)
out <- file.path(root, "results/p6/supplementary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(row, file.path(out, "scz_correction_off_paired_phbc.tsv"), sep = "\t")
cat("P6_SCZ_CORRECTION_OFF_COMPLETE\n")

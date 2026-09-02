#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages(library(data.table))
records <- readRDS(file.path(root, "results/p5/neurological_phbc/formal/neurological_phbc_records.rds"))
n_block <- 200L
left_key <- "EPILEPSY_EHR_META_2OF3"
right_key <- "EPILEPSY_EHR_META_3OF3"
left <- records[[left_key]]
right <- records[[right_key]]
stopifnot(identical(left$auxiliaries, right$auxiliaries),
          length(left$block_uncorrected) == n_block,
          length(right$block_uncorrected) == n_block)
delta <- left$phbc - right$phbc
delta_jk <- left$corrected_weight^2 * left$block_uncorrected -
  right$corrected_weight^2 * right$block_uncorrected
pseudo <- n_block * delta - (n_block - 1L) * delta_jk
se <- sd(pseudo) / sqrt(n_block)
row <- data.table(
  comparison = "EPILEPSY_EHR_META_2OF3_minus_3OF3",
  analysis_tier = "meta_definition_sensitivity",
  left_target = left_key,
  right_target = right_key,
  left_phbc = left$phbc,
  right_phbc = right$phbc,
  delta_pp = 100 * delta,
  paired_jackknife_se_pp = 100 * se,
  ci95_low_pp = 100 * (delta - 1.96 * se),
  ci95_high_pp = 100 * (delta + 1.96 * se),
  z = delta / se,
  p_two_sided = 2 * pnorm(-abs(delta / se)),
  n_block = n_block,
  left_corrected_weight = left$corrected_weight,
  right_corrected_weight = right$corrected_weight,
  method = "200-block paired delete-one-block jackknife; meta-definition sensitivity"
)
out <- file.path(root, "results/p6/supplementary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(row, file.path(out, "epilepsy_ehr_meta_2of3_minus_3of3_paired_phbc.tsv"), sep = "\t")
cat("P6_EPILEPSY_META_COMPOSITION_PAIR_COMPLETE\n")

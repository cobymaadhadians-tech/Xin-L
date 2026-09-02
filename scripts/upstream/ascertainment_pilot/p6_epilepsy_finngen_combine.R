#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({library(data.table); library(pleioh2g)})
target <- "EPILEPSY_FINNGEN_R13"
aux <- c("NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641",
         "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387")
traits <- c(target, aux)
n_block <- 200L
sample_rep <- 1000L
files <- sort(list.files(file.path(root, "results/p6/epilepsy_finngen/cache"),
                         pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
if (length(files) != 20L) stop("expected 20 FinnGen cache chunks")
chunks <- lapply(files, readRDS)
ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) stop("incomplete cache")
rgarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
h2array <- matrix(NA_real_, n_block, length(traits), dimnames = list(paste0("block_", seq_len(n_block)), traits))
gcovarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
for (x in chunks) {
  rgarray[, , x$block_ids] <- x$rgarray
  h2array[x$block_ids, ] <- x$h2array
  gcovarray[, , x$block_ids] <- x$gcovarray
}
full <- readRDS(file.path(root, "results/p6/epilepsy_finngen/full_matrix/ldsc_alltraits_observed.rds"))
full_rg <- full$rg[traits, traits, drop = FALSE]
full_h2 <- full$h2[, traits, drop = FALSE]
aux_rg <- full_rg[aux, aux, drop = FALSE]
eig <- eigen(aux_rg, symmetric = TRUE, only.values = TRUE)$values
target_aux_max_abs_rg <- max(abs(full_rg[target, aux]))
if (min(eig) <= 0 || target_aux_max_abs_rg >= sqrt(0.5)) stop("FinnGen neurological panel failed matrix or target-auxiliary gate")
block_uncorrected <- vapply(seq_len(n_block), function(b) {
  r <- rgarray[1L, -1L, b]
  R <- rgarray[-1L, -1L, b]
  as.numeric(crossprod(r, solve(R, r)))
}, numeric(1L))
set.seed(20260831)
result <- pleiotropyh2_cor_computing_single(1L, traits, full_h2, h2array, full_rg, rgarray, sample_rep)
phbc <- as.numeric(result$percentage_h2pleio_corr)
weight <- as.numeric(result$corrected_weight)
official_se <- as.numeric(result$percentage_h2pleio_corr_se)
official_z <- as.numeric(result$percentage_h2pleio_corr_Z)
if (!all(is.finite(c(phbc, weight, official_se, official_z))) || weight < 0.5 || official_se > 0.5) stop("FinnGen PHBC correction stability failure")

ilae <- readRDS(file.path(root, "results/p5/neurological_phbc/formal/neurological_phbc_records.rds"))[["EPILEPSY_ILAE2023_EUR"]]
delta <- phbc - ilae$phbc
delta_jk <- weight^2 * block_uncorrected - ilae$corrected_weight^2 * ilae$block_uncorrected
pseudo <- n_block * delta - (n_block - 1L) * delta_jk
se <- sd(pseudo) / sqrt(n_block)
out <- file.path(root, "results/p6/epilepsy_finngen")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(target = target, auxiliaries = paste(aux, collapse = ","), n_aux = length(aux),
                  phbc = phbc, phbc_pp = 100 * phbc, official_se = official_se,
                  official_z = official_z, corrected_weight = weight,
                  block_jackknife_se = sd(n_block * phbc - (n_block - 1L) * weight^2 * block_uncorrected) / sqrt(n_block),
                  target_aux_max_abs_rg = target_aux_max_abs_rg,
                  aux_min_eigenvalue = min(eig), aux_condition_number = max(eig) / min(eig), gate = "PASS"),
       file.path(out, "finngen_neurological_panel_phbc.tsv"), sep = "\t")
fwrite(data.table(comparison = "EPILEPSY_FINNGEN_minus_ILAE", analysis_tier = "neurological_sensitivity_direct_registry",
                  left_target = target, right_target = "EPILEPSY_ILAE2023_EUR",
                  left_phbc = phbc, right_phbc = ilae$phbc, delta_pp = 100 * delta,
                  paired_jackknife_se_pp = 100 * se,
                  ci95_low_pp = 100 * (delta - 1.96 * se), ci95_high_pp = 100 * (delta + 1.96 * se),
                  z = delta / se, p_two_sided = 2 * pnorm(-abs(delta / se)), n_block = n_block,
                  left_corrected_weight = weight, right_corrected_weight = ilae$corrected_weight),
       file.path(out, "finngen_minus_ilae_paired_phbc.tsv"), sep = "\t")
saveRDS(list(target = target, auxiliaries = aux, full_rg = full_rg, full_h2 = full_h2,
             rgarray = rgarray, h2array = h2array, gcovarray = gcovarray,
             phbc = phbc, corrected_weight = weight, block_uncorrected = block_uncorrected,
             result = result), file.path(out, "finngen_neurological_panel_phbc_record.rds"))
cat("P6_EPILEPSY_FINNGEN_COMBINE_COMPLETE\n")

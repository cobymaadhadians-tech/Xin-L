#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
n_block <- 200L
sample_rep <- 1000L
targets <- c(
  "AD_GCST90704646_MAIN", "AD_GCST90704647_NOPROXY", "AD_GCST90704648_NOBIOBANK",
  "EPILEPSY_ILAE2023_EUR", "EPILEPSY_EHR_META_2OF3", "EPILEPSY_EHR_META_3OF3"
)
aux <- c("NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641", "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387")
out_dir <- file.path(root, "results/p5/neurological_phbc/formal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

uncorrected_blocks <- function(rgarray) {
  vapply(seq_len(dim(rgarray)[3L]), function(b) {
    r <- rgarray[1L, -1L, b]
    R <- rgarray[-1L, -1L, b]
    as.numeric(t(r) %*% solve(R, r))
  }, numeric(1L))
}

records <- setNames(vector("list", length(targets)), targets)
summary_rows <- list()
for (target in targets) {
  traits <- c(target, aux)
  chunk_dir <- file.path(root, "results/p5/neurological_phbc/cache", target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, function(x) x$block_ids), use.names = FALSE)
  if (!identical(sort(ids), seq_len(n_block)) || anyDuplicated(ids)) stop("incomplete cache: ", target)
  rgarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  h2array <- matrix(NA_real_, n_block, length(traits), dimnames = list(paste0("block_", seq_len(n_block)), traits))
  gcovarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  for (x in chunks) {
    if (!identical(x$traits, traits)) stop("trait mismatch: ", target)
    rgarray[, , x$block_ids] <- x$rgarray
    h2array[x$block_ids, ] <- x$h2array
    gcovarray[, , x$block_ids] <- x$gcovarray
  }
  full <- readRDS(file.path(root, "results/p5/neurological_phbc/full_matrix", target, "ldsc_alltraits_observed.rds"))
  full_h2 <- full$h2[, traits, drop = FALSE]
  full_rg <- full$rg[traits, traits, drop = FALSE]
  block_uncorr <- uncorrected_blocks(rgarray)
  set.seed(20260830)
  result <- tryCatch(
    pleiotropyh2_cor_computing_single(1L, traits, full_h2, h2array, full_rg, rgarray, sample_rep),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    r <- full_rg[1L, -1L]
    phbc <- as.numeric(t(r) %*% solve(full_rg[-1L, -1L], r))
    weight <- 1
    official_se <- NA_real_
    official_z <- NA_real_
    method <- paste0("uncorrected quadratic form; official correction failed: ", conditionMessage(result))
    result_saved <- NULL
  } else {
    phbc <- as.numeric(result$percentage_h2pleio_corr)
    weight <- as.numeric(result$corrected_weight)
    official_se <- as.numeric(result$percentage_h2pleio_corr_se)
    official_z <- as.numeric(result$percentage_h2pleio_corr_Z)
    method <- "official direct correction on frozen five-trait neurological panel"
    result_saved <- result
  }
  pseudo <- n_block * phbc - (n_block - 1L) * weight^2 * block_uncorr
  block_se <- sd(pseudo) / sqrt(n_block)
  records[[target]] <- list(target = target, traits = traits, auxiliaries = aux, n_block = n_block,
    full_h2 = full_h2, full_rg = full_rg, h2array = h2array, rgarray = rgarray, gcovarray = gcovarray,
    phbc = phbc, corrected_weight = weight, block_uncorrected = block_uncorr, block_se = block_se,
    official_se = official_se, official_z = official_z, method = method, result = result_saved)
  summary_rows[[length(summary_rows) + 1L]] <- data.table(target = target, n_aux = length(aux), auxiliaries = paste(aux, collapse = ","),
    phbc = phbc, phbc_pp = 100 * phbc, official_se = official_se, official_z = official_z,
    block_jackknife_se = block_se, corrected_weight = weight, method = method)
}

pairs <- data.table(
  comparison = c("AD_NOPROXY_minus_MAIN", "AD_NOBIOBANK_minus_MAIN", "EPILEPSY_EHR2_minus_ILAE", "EPILEPSY_EHR3_minus_ILAE"),
  family = c("primary", "sensitivity", "primary", "sensitivity"),
  left = c("AD_GCST90704647_NOPROXY", "AD_GCST90704648_NOBIOBANK", "EPILEPSY_EHR_META_2OF3", "EPILEPSY_EHR_META_3OF3"),
  right = c("AD_GCST90704646_MAIN", "AD_GCST90704646_MAIN", "EPILEPSY_ILAE2023_EUR", "EPILEPSY_ILAE2023_EUR")
)
delta_rows <- list()
for (i in seq_len(nrow(pairs))) {
  left <- records[[pairs$left[[i]]]]
  right <- records[[pairs$right[[i]]]]
  delta <- left$phbc - right$phbc
  delta_jk <- left$corrected_weight^2 * left$block_uncorrected - right$corrected_weight^2 * right$block_uncorrected
  pseudo <- n_block * delta - (n_block - 1L) * delta_jk
  se <- sd(pseudo) / sqrt(n_block)
  delta_rows[[i]] <- data.table(comparison = pairs$comparison[[i]], family = pairs$family[[i]], left_target = pairs$left[[i]], right_target = pairs$right[[i]],
    delta_pp = 100 * delta, paired_jackknife_se_pp = 100 * se, ci95_low_pp = 100 * (delta - 1.96 * se),
    ci95_high_pp = 100 * (delta + 1.96 * se), z = delta / se, p_two_sided = 2 * pnorm(-abs(delta / se)),
    left_method = left$method, right_method = right$method)
}
delta <- rbindlist(delta_rows)
delta[, q_bh_within_family := p.adjust(p_two_sided, method = "BH"), by = family]
delta[, survives_bh_0_05 := q_bh_within_family < 0.05]

fwrite(rbindlist(summary_rows), file.path(out_dir, "neurological_phbc_summary.tsv"), sep = "\t")
fwrite(delta, file.path(out_dir, "neurological_paired_delta_phbc.tsv"), sep = "\t")
saveRDS(records, file.path(out_dir, "neurological_phbc_records.rds"))
cat("P5_NEURO_PHBC_COMBINE_COMPLETE\n")

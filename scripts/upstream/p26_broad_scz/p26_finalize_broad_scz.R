#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})

n_block <- 200L
sample_rep <- 1000L
seed <- 20260831L
broad <- "SCZ_FINNGEN_R13_EXMORE"
comparators <- c("SCZ_PGC2022", "SCZ_FINNGEN_R13")
aux <- c(
  "MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD", "ASD",
  "OCD_2025", "AN", "AUD", "CUD_2023_EUR"
)
domains <- list(
  mood_psychotic = c("MDD_CLIN_PGC2025", "BD_PGC2021"),
  internalizing = c("PTSD", "ANX"),
  neurodevelopmental = c("ADHD", "ASD"),
  compulsive_eating = c("OCD_2025", "AN"),
  substance = c("AUD", "CUD_2023_EUR")
)
out_dir <- file.path(root, "results/p26_broad_scz/final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combine_chunks <- function(target, traits, cache_root) {
  chunk_dir <- file.path(cache_root, target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop("expected 20 chunks for ", target, "; found ", length(files))
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) stop("invalid block coverage for ", target)
  rgarray <- array(NA_real_, c(length(traits), length(traits), n_block),
                   dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  h2array <- matrix(NA_real_, n_block, length(traits),
                    dimnames = list(paste0("block_", seq_len(n_block)), traits))
  gcovarray <- array(NA_real_, c(length(traits), length(traits), n_block),
                     dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  for (x in chunks) {
    if (!identical(x$target, target) || !identical(x$traits, traits)) stop("trait order mismatch: ", target)
    idx <- x$block_ids
    rgarray[, , idx] <- x$rgarray
    h2array[idx, ] <- x$h2array
    gcovarray[, , idx] <- x$gcovarray
  }
  if (!all(is.finite(rgarray)) || !all(is.finite(h2array)) || !all(is.finite(gcovarray))) {
    stop("non-finite cache values: ", target)
  }
  list(rgarray = rgarray, h2array = h2array, gcovarray = gcovarray)
}

load_base <- function(target) {
  traits <- c(target, aux)
  if (target == broad) {
    matrix_root <- file.path(root, "results/p26_broad_scz/full_matrix")
    cache_root <- file.path(root, "results/p26_broad_scz/phbc_cache")
  } else {
    matrix_root <- file.path(root, "results/p25/full_matrix")
    cache_root <- file.path(root, "results/p25/phbc_cache")
  }
  full <- readRDS(file.path(matrix_root, target, "ldsc_alltraits_observed.rds"))
  cache <- combine_chunks(target, traits, cache_root)
  list(
    target = target, traits = traits, auxiliaries = aux,
    full_h2 = full$h2[, traits, drop = FALSE],
    full_rg = full$rg[traits, traits, drop = FALSE],
    h2array = cache$h2array, rgarray = cache$rgarray, gcovarray = cache$gcovarray
  )
}

subset_input <- function(x, auxiliaries) {
  traits <- c(x$target, auxiliaries)
  idx <- match(traits, x$traits)
  if (anyNA(idx)) stop("missing trait in subset: ", x$target)
  list(
    target = x$target, traits = traits, auxiliaries = auxiliaries,
    full_h2 = x$full_h2[, idx, drop = FALSE],
    full_rg = x$full_rg[idx, idx, drop = FALSE],
    h2array = x$h2array[, idx, drop = FALSE],
    rgarray = x$rgarray[idx, idx, , drop = FALSE],
    gcovarray = x$gcovarray[idx, idx, , drop = FALSE]
  )
}

uncorrected_point <- function(x) {
  r <- x$full_rg[1L, -1L]
  r_aux <- x$full_rg[-1L, -1L, drop = FALSE]
  as.numeric(crossprod(r, solve(r_aux, r)))
}

uncorrected_blocks <- function(x) {
  vapply(seq_len(n_block), function(b) {
    r <- x$rgarray[1L, -1L, b]
    r_aux <- x$rgarray[-1L, -1L, b, drop = FALSE]
    dim(r_aux) <- c(length(r), length(r))
    as.numeric(crossprod(r, solve(r_aux, r)))
  }, numeric(1L))
}

fit_full <- function(x) {
  set.seed(seed)
  result <- pleiotropyh2_cor_computing_single(
    1L, x$traits, x$full_h2, x$h2array, x$full_rg, x$rgarray, sample_rep
  )
  weight <- as.numeric(result$corrected_weight)
  point_uncorrected <- uncorrected_point(x)
  blocks_uncorrected <- uncorrected_blocks(x)
  list(
    target = x$target,
    auxiliaries = x$auxiliaries,
    phbc = as.numeric(result$percentage_h2pleio_corr),
    official_se = as.numeric(result$percentage_h2pleio_corr_se),
    official_z = as.numeric(result$percentage_h2pleio_corr_Z),
    corrected_weight = weight,
    uncorrected = point_uncorrected,
    block_uncorrected = blocks_uncorrected,
    corrected_blocks = weight^2 * blocks_uncorrected,
    selected_auxD = as.character(result$selected_auxD)
  )
}

paired_test <- function(point, blocks) {
  pseudo <- n_block * point - (n_block - 1L) * blocks
  se <- sd(pseudo) / sqrt(n_block)
  z <- point / se
  list(se = se, low = point - 1.96 * se, high = point + 1.96 * se,
       z = z, p = 2 * pnorm(-abs(z)))
}

base <- setNames(lapply(c(broad, comparators), load_base), c(broad, comparators))
records_path <- file.path(root, "results/p11_scz_neff_audit/corrected_p25_final/p25_phbc_records.rds")
stored <- readRDS(records_path)$full
full_records <- list()
full_records[[broad]] <- fit_full(base[[broad]])
for (target in comparators) {
  x <- stored[[target]]
  if (is.null(x)) stop("missing stored comparator record: ", target)
  full_records[[target]] <- list(
    target = target, auxiliaries = x$auxiliaries, phbc = x$phbc,
    official_se = x$official_se, official_z = x$official_z,
    corrected_weight = x$corrected_weight,
    uncorrected = x$uncorrected_quadratic_form,
    block_uncorrected = x$block_uncorrected,
    corrected_blocks = x$corrected_blocks,
    selected_auxD = x$selected_auxD
  )
}

full_summary <- rbindlist(lapply(full_records, function(x) data.table(
  target = x$target, auxiliaries = paste(x$auxiliaries, collapse = ","),
  n_aux = length(x$auxiliaries), phbc = x$phbc, phbc_pp = 100 * x$phbc,
  official_se = x$official_se, official_z = x$official_z,
  corrected_weight = x$corrected_weight,
  uncorrected_quadratic_form = x$uncorrected,
  selected_auxD = paste(x$selected_auxD, collapse = ",")
)))

comparison_rows <- lapply(comparators, function(right) {
  left_rec <- full_records[[broad]]
  right_rec <- full_records[[right]]
  point <- left_rec$phbc - right_rec$phbc
  tst <- paired_test(point, left_rec$corrected_blocks - right_rec$corrected_blocks)
  raw_point <- left_rec$uncorrected - right_rec$uncorrected
  raw_tst <- paired_test(raw_point, left_rec$block_uncorrected - right_rec$block_uncorrected)
  data.table(
    comparison = paste(broad, "minus", right), left_target = broad, right_target = right,
    left_phbc_pp = 100 * left_rec$phbc, right_phbc_pp = 100 * right_rec$phbc,
    delta_pp = 100 * point, paired_se_pp = 100 * tst$se,
    ci95_low_pp = 100 * tst$low, ci95_high_pp = 100 * tst$high,
    z = tst$z, p_two_sided = tst$p,
    correction_off_delta_pp = 100 * raw_point,
    correction_off_se_pp = 100 * raw_tst$se,
    correction_off_ci95_low_pp = 100 * raw_tst$low,
    correction_off_ci95_high_pp = 100 * raw_tst$high,
    correction_off_p_two_sided = raw_tst$p
  )
})
paired_total <- rbindlist(comparison_rows)
paired_total[, q_bh_2 := p.adjust(p_two_sided, method = "BH")]

category_rows <- list()
category_records <- list()
for (target in c(broad, comparators)) {
  full_rec <- full_records[[target]]
  for (domain in names(domains)) {
    remaining <- setdiff(aux, domains[[domain]])
    reduced <- subset_input(base[[target]], remaining)
    reduced_uncorrected <- uncorrected_point(reduced)
    reduced_blocks <- uncorrected_blocks(reduced)
    scale <- full_rec$corrected_weight^2
    reduction <- scale * (full_rec$uncorrected - reduced_uncorrected)
    reduction_blocks <- scale * (full_rec$block_uncorrected - reduced_blocks)
    key <- paste(target, domain, sep = "__")
    category_records[[key]] <- list(point = reduction, blocks = reduction_blocks)
    category_rows[[length(category_rows) + 1L]] <- data.table(
      target = target, domain = domain,
      removed_auxiliaries = paste(domains[[domain]], collapse = ","),
      reduction_pp = 100 * reduction, full_panel_scale = scale
    )
  }
}
category_estimates <- rbindlist(category_rows)

category_delta_rows <- list()
for (right in comparators) {
  for (domain in names(domains)) {
    left_rec <- category_records[[paste(broad, domain, sep = "__")]]
    right_rec <- category_records[[paste(right, domain, sep = "__")]]
    point <- left_rec$point - right_rec$point
    tst <- paired_test(point, left_rec$blocks - right_rec$blocks)
    category_delta_rows[[length(category_delta_rows) + 1L]] <- data.table(
      comparison = paste(broad, "minus", right), domain = domain,
      broad_reduction_pp = 100 * left_rec$point,
      comparator_reduction_pp = 100 * right_rec$point,
      difference_in_reduction_pp = 100 * point,
      paired_se_pp = 100 * tst$se,
      ci95_low_pp = 100 * tst$low, ci95_high_pp = 100 * tst$high,
      z = tst$z, p_two_sided = tst$p
    )
  }
}
category_delta <- rbindlist(category_delta_rows)
category_delta[, q_bh_10 := p.adjust(p_two_sided, method = "BH")]

profile_rows <- list()
profile_tests <- list()
for (right in comparators) {
  point <- base[[broad]]$full_rg[broad, aux] - base[[right]]$full_rg[right, aux]
  block_diff <- t(vapply(seq_len(n_block), function(b) {
    base[[broad]]$rgarray[broad, aux, b] - base[[right]]$rgarray[right, aux, b]
  }, numeric(length(aux))))
  pseudo <- sweep(block_diff, 2L, n_block * point, function(x, y) y - (n_block - 1L) * x)
  cov_est <- cov(pseudo) / n_block
  se_each <- sqrt(diag(cov_est))
  z_each <- point / se_each
  comparison <- paste(broad, "minus", right)
  profile_rows[[comparison]] <- data.table(
    comparison = comparison, auxiliary = aux, delta_rg = as.numeric(point),
    paired_se = se_each, z = z_each, p_two_sided = 2 * pnorm(-abs(z_each))
  )
  w <- rep(1 / length(aux), length(aux))
  mean_delta <- mean(point)
  mean_se <- sqrt(as.numeric(t(w) %*% cov_est %*% w))
  omnibus <- as.numeric(t(point) %*% solve(cov_est, point))
  profile_tests[[comparison]] <- data.table(
    comparison = comparison, mean_delta_rg = mean_delta,
    mean_paired_se = mean_se, mean_z = mean_delta / mean_se,
    mean_p_two_sided = 2 * pnorm(-abs(mean_delta / mean_se)),
    omnibus_chisq = omnibus, omnibus_df = length(aux),
    omnibus_p = pchisq(omnibus, df = length(aux), lower.tail = FALSE),
    covariance_condition_number = kappa(cov_est)
  )
}
profile <- rbindlist(profile_rows)
profile_test <- rbindlist(profile_tests)

fwrite(full_summary, file.path(out_dir, "p26_scz_full_panel_phbc.tsv"), sep = "\t")
fwrite(paired_total, file.path(out_dir, "p26_scz_paired_total_phbc.tsv"), sep = "\t")
fwrite(category_estimates, file.path(out_dir, "p26_scz_category_reductions_official.tsv"), sep = "\t")
fwrite(category_delta, file.path(out_dir, "p26_scz_category_paired_did_official.tsv"), sep = "\t")
fwrite(profile, file.path(out_dir, "p26_scz_rg_profile_paired.tsv"), sep = "\t")
fwrite(profile_test, file.path(out_dir, "p26_scz_rg_profile_tests.tsv"), sep = "\t")
saveRDS(list(full = full_records, category = category_records),
        file.path(out_dir, "p26_scz_phbc_records.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines("PASS", file.path(out_dir, "STATUS.txt"))
cat("P26_FINALIZE_COMPLETE\n")

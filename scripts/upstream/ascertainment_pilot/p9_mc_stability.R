#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(mvtnorm)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

n_block <- 200L
targets <- c(
  "MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025",
  "SCZ_PGC2022", "SCZ_FINNGEN_R13"
)
settings <- expand.grid(
  sample_rep = c(1000L, 10000L),
  seed = c(123L, 20260831L, 8675309L),
  KEEP.OUT.ATTRS = FALSE
)
out_dir <- file.path(root, "results/p9_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lower_tri_vector <- function(mat) mat[which(lower.tri(mat, diag = FALSE), arr.ind = TRUE)]

load_target <- function(target) {
  traits <- c(target, p25_panels[[target]])
  full <- readRDS(file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds"))
  files <- sort(list.files(
    file.path(root, "results/p25/phbc_cache", target),
    pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE
  ))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target)
  chunks <- lapply(files, readRDS)
  block_ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(block_ids) || !identical(sort(block_ids), seq_len(n_block))) {
    stop("incomplete or duplicated blocks for ", target)
  }
  rgarray <- array(
    NA_real_, c(length(traits), length(traits), n_block),
    dimnames = list(traits, traits, paste0("block_", seq_len(n_block)))
  )
  for (chunk in chunks) {
    if (!identical(chunk$traits, traits)) stop("trait-order mismatch for ", target)
    rgarray[, , chunk$block_ids] <- chunk$rgarray
  }
  values <- vapply(seq_len(n_block), function(b) lower_tri_vector(rgarray[, , b]), numeric(length(traits) * (length(traits) - 1L) / 2L))
  list(
    target = target,
    traits = traits,
    point_rg = full$rg[traits, traits],
    covariance = stats::cov(t(values)) * n_block
  )
}

draw_noise <- function(point_rg, covariance, n_draw, seed) {
  set.seed(seed)
  n <- nrow(point_rg)
  out <- vector("list", n_draw)
  attempts <- integer(n_draw)
  for (sample_id in seq_len(n_draw)) {
    for (iter in seq_len(202L)) {
      lower <- matrix(0, nrow = n, ncol = n)
      lower[lower.tri(lower, diag = FALSE)] <- mvtnorm::rmvnorm(1L, sigma = covariance)
      noise <- lower + t(lower)
      noisy_rg <- point_rg + noise
      if (all(eigen(noisy_rg, symmetric = TRUE, only.values = TRUE)$values > 0)) {
        out[[sample_id]] <- noise
        attempts[[sample_id]] <- iter
        break
      }
    }
    if (is.null(out[[sample_id]])) stop("failed to draw a positive-definite matrix")
  }
  list(noises = out, attempts = attempts)
}

sample_phbc <- function(point_rg, target_index, noises, ratio) {
  target_rg <- point_rg[target_index, -target_index]
  vapply(noises, function(noise) {
    noisy_target <- ratio * target_rg + noise[target_index, -target_index]
    noisy_aux <- point_rg[-target_index, -target_index, drop = FALSE] +
      noise[-target_index, -target_index, drop = FALSE]
    as.numeric(crossprod(noisy_target, solve(noisy_aux, noisy_target)))
  }, numeric(1L))
}

estimate_weight <- function(x, n_draw, seed) {
  target_index <- match(x$target, x$traits)
  truth <- as.numeric(pleioh2g:::Cal_cor_test_single(x$point_rg, target_index))
  draw <- draw_noise(x$point_rg, x$covariance, n_draw, seed)
  initial <- sample_phbc(x$point_rg, target_index, draw$noises, 1)
  initial_mean <- mean(initial)
  initial_sd <- stats::sd(initial)
  lower <- 0
  upper <- 1
  ratio <- 1
  current_mean <- initial_mean
  current_sd <- initial_sd
  rounds <- 0L
  while ((upper - lower) * max(abs(x$point_rg[target_index, -target_index])) > 1e-6) {
    rounds <- rounds + 1L
    if (current_mean < truth) lower <- ratio else upper <- ratio
    ratio <- (lower + upper) / 2
    values <- sample_phbc(x$point_rg, target_index, draw$noises, ratio)
    current_mean <- mean(values)
    current_sd <- stats::sd(values)
    if (ratio < 0.5) stop("correction weight fell below 0.5")
    if (abs(current_mean - truth) / truth < 0.05) break
  }
  data.frame(
    target = x$target,
    sample_rep = n_draw,
    seed = seed,
    uncorrected_phbc = truth,
    corrected_weight = ratio,
    corrected_phbc = ratio^2 * truth,
    calibrated_relative_error = abs(current_mean - truth) / truth,
    mc_sd_ratio = current_sd / initial_sd,
    binary_search_rounds = rounds,
    max_pd_attempts = max(draw$attempts),
    stringsAsFactors = FALSE
  )
}

inputs <- setNames(lapply(targets, load_target), targets)
rows <- list()
for (target in targets) {
  for (i in seq_len(nrow(settings))) {
    rows[[length(rows) + 1L]] <- estimate_weight(
      inputs[[target]], settings$sample_rep[[i]], settings$seed[[i]]
    )
  }
}
target_results <- do.call(rbind, rows)

pairs <- data.frame(
  comparison = c(
    "MDD_EHR_minus_QUEST", "MDD_FINNGEN_minus_QUEST",
    "SCZ_FINNGEN_minus_PGC"
  ),
  left = c("MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "SCZ_FINNGEN_R13"),
  right = c("MDD_QUEST_PGC2025", "MDD_QUEST_PGC2025", "SCZ_PGC2022"),
  stringsAsFactors = FALSE
)

contrast_rows <- list()
for (i in seq_len(nrow(settings))) {
  subset <- target_results[
    target_results$sample_rep == settings$sample_rep[[i]] &
      target_results$seed == settings$seed[[i]],
  ]
  for (j in seq_len(nrow(pairs))) {
    left <- subset[subset$target == pairs$left[[j]], ]
    right <- subset[subset$target == pairs$right[[j]], ]
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      comparison = pairs$comparison[[j]],
      sample_rep = settings$sample_rep[[i]],
      seed = settings$seed[[i]],
      left_corrected_weight = left$corrected_weight,
      right_corrected_weight = right$corrected_weight,
      delta_fraction = left$corrected_phbc - right$corrected_phbc,
      delta_pp = 100 * (left$corrected_phbc - right$corrected_phbc),
      stringsAsFactors = FALSE
    )
  }
}
contrast_results <- do.call(rbind, contrast_rows)

write.table(
  target_results, file.path(out_dir, "mc_stability_target_phbc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  contrast_results, file.path(out_dir, "mc_stability_headline_contrasts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("P9_MC_STABILITY_COMPLETE\n")

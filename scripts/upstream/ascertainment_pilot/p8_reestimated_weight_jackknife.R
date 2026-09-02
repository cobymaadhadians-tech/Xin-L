#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("usage: p8_reestimated_weight_jackknife.R TARGET START_BLOCK END_BLOCK [LABEL]")
}

target <- args[[1L]]
start_block <- as.integer(args[[2L]])
end_block <- as.integer(args[[3L]])
label <- if (length(args) >= 4L) args[[4L]] else sprintf("%03d_%03d", start_block, end_block)

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(mvtnorm)
  library(pleioh2g)
})

sample_rep <- 1000L
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))
out_dir <- file.path(root, "results/p8/p25_reestimated_weight_jackknife", target)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!target %in% names(p25_panels)) stop("target is absent from the P25 panel configuration")
traits <- c(target, p25_panels[[target]])
full <- readRDS(file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds"))
chunk_dir <- file.path(root, "results/p25/phbc_cache", target)
chunk_files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
if (length(chunk_files) != 20L) stop("expected 20 P25 cache chunks for ", target)
chunks <- lapply(chunk_files, readRDS)
block_ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
if (anyDuplicated(block_ids) || !identical(sort(block_ids), seq_len(200L))) {
  stop("incomplete or duplicated P25 blocks for ", target)
}
rgarray <- array(
  NA_real_, c(length(traits), length(traits), 200L),
  dimnames = list(traits, traits, paste0("block_", seq_len(200L)))
)
for (chunk in chunks) {
  if (!identical(chunk$traits, traits)) stop("P25 cache trait-order mismatch for ", target)
  rgarray[, , chunk$block_ids] <- chunk$rgarray
}
records <- readRDS(file.path(root, "results/p25/final/p25_phbc_records.rds"))
record <- records$full[[target]]
x <- list(
  target = target,
  auxiliaries = p25_panels[[target]],
  n_block = 200L,
  full_rg = full$rg[traits, traits, drop = FALSE],
  rgarray = rgarray,
  result = record$result
)
n_block <- as.integer(x$n_block)
if (start_block < 1L || end_block > n_block || start_block > end_block) {
  stop("invalid block range")
}

lower_tri_vector <- function(mat) {
  mat[which(lower.tri(mat, diag = FALSE), arr.ind = TRUE)]
}

jackknife_noise_covariance <- function(rgarray) {
  n <- dim(rgarray)[1L]
  b <- dim(rgarray)[3L]
  values <- matrix(NA_real_, nrow = n * (n - 1L) / 2L, ncol = b)
  for (i in seq_len(b)) values[, i] <- lower_tri_vector(rgarray[, , i])
  stats::cov(t(values)) * b
}

draw_accepted_noise <- function(point_rg, covariance, n_draw) {
  n <- nrow(point_rg)
  out <- vector("list", n_draw)
  set.seed(123)
  for (sample_id in seq_len(n_draw)) {
    accepted <- FALSE
    for (iter in seq_len(202L)) {
      lower <- matrix(0, nrow = n, ncol = n)
      lower[lower.tri(lower, diag = FALSE)] <- mvtnorm::rmvnorm(1L, sigma = covariance)
      noise <- lower + t(lower)
      noisy_rg <- point_rg + noise
      if (all(eigen(noisy_rg, symmetric = TRUE, only.values = TRUE)$values > 0) &&
          all(diag(noisy_rg) > 0)) {
        out[[sample_id]] <- noise
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("failed to draw a positive-definite matrix for sample ", sample_id)
  }
  out
}

sample_phbc <- function(point_rg, target_index, noises, ratio) {
  target_rg <- point_rg[target_index, -target_index]
  vapply(noises, function(noise) {
    noisy_target <- ratio * target_rg + noise[target_index, -target_index]
    noisy_aux <- point_rg[-target_index, -target_index, drop = FALSE] +
      noise[-target_index, -target_index, drop = FALSE]
    as.numeric(t(noisy_target) %*% solve(noisy_aux, noisy_target))
  }, numeric(1L))
}

estimate_weight <- function(point_rg, target_index, covariance, n_draw = 1000L) {
  truth <- as.numeric(
    pleioh2g:::Cal_cor_test_single(point_rg, target_index)
  )
  noises <- draw_accepted_noise(point_rg, covariance, n_draw)
  initial <- sample_phbc(point_rg, target_index, noises, 1)
  initial_mean <- mean(initial)
  initial_sd <- stats::sd(initial)

  lower <- 0
  upper <- 1
  ratio <- 1
  current_mean <- initial_mean
  current_sd <- initial_sd
  round <- 0L

  while ((upper - lower) * max(abs(point_rg[target_index, -target_index])) > 1e-6) {
    round <- round + 1L
    if (current_mean < truth) lower <- ratio else upper <- ratio
    ratio <- (lower + upper) / 2
    values <- sample_phbc(point_rg, target_index, noises, ratio)
    current_mean <- mean(values)
    current_sd <- stats::sd(values)
    if (ratio < 0.5) stop("correction weight fell below 0.5")
    if (abs(current_mean - truth) / truth < 0.05) break
  }

  list(
    weight = ratio,
    truth = truth,
    initial_mean = initial_mean,
    calibrated_mean = current_mean,
    initial_sd = initial_sd,
    calibrated_sd = current_sd,
    relative_error = abs(current_mean - truth) / truth,
    rounds = round
  )
}

target_index <- match(target, dimnames(x$rgarray)[[1L]])
if (is.na(target_index)) stop("target is absent from rgarray dimnames")
covariance <- jackknife_noise_covariance(x$rgarray)

full_check <- estimate_weight(x$full_rg, target_index, covariance, sample_rep)
stored_weight <- as.numeric(x$result$corrected_weight)
if (!isTRUE(all.equal(full_check$weight, stored_weight, tolerance = 1e-12))) {
  stop("full-data correction weight mismatch: recalculated=", full_check$weight,
       ", stored=", stored_weight)
}

rows <- vector("list", end_block - start_block + 1L)
for (j in seq_along(rows)) {
  block <- start_block + j - 1L
  point_rg <- x$rgarray[, , block]
  estimate <- estimate_weight(point_rg, target_index, covariance, sample_rep)
  uncorrected <- as.numeric(
    pleioh2g:::Cal_cor_test_single(point_rg, target_index)
  )
  rows[[j]] <- data.frame(
    target = target,
    block = block,
    reestimated_weight = estimate$weight,
    fixed_weight = stored_weight,
    reestimated_phbc = estimate$weight^2 * uncorrected,
    fixed_weight_phbc = stored_weight^2 * uncorrected,
    uncorrected_phbc = uncorrected,
    calibrated_relative_error = estimate$relative_error,
    binary_search_rounds = estimate$rounds,
    sample_rep = sample_rep,
    covariance_blocks = n_block,
    stringsAsFactors = FALSE
  )
  cat(target, "block", block, "weight", estimate$weight, "\n")
}

out <- do.call(rbind, rows)
write.table(
  out,
  file.path(out_dir, paste0("blocks_", label, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(
  list(
    target = target,
    block_range = c(start_block, end_block),
    sample_rep = sample_rep,
    covariance_blocks = n_block,
    full_weight_check = full_check,
    stored_full_weight = stored_weight,
    results = out
  ),
  file.path(out_dir, paste0("blocks_", label, ".rds"))
)

cat("P8_REESTIMATED_WEIGHT_CHUNK_COMPLETE\n")

#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

n_block <- 200L
sample_rep <- 1000L
seed <- 20260831L
out_dir <- file.path(root, "results/p11_scz_neff_audit/corrected_leave_one_auxiliary_out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

headline_pairs <- data.table(
  comparison = c(
    "SCZ_FINNGEN_minus_PGC",
    "MDD_EHR_minus_QUEST",
    "MDD_FINNGEN_minus_QUEST"
  ),
  disease = c("SCZ", "MDD", "MDD"),
  left = c("SCZ_FINNGEN_R13", "MDD_EHR_PGC2025", "MDD_FINNGEN_R13"),
  right = c("SCZ_PGC2022", "MDD_QUEST_PGC2025", "MDD_QUEST_PGC2025")
)

required_targets <- unique(c(headline_pairs$left, headline_pairs$right))

combine_chunks <- function(target, traits) {
  chunk_dir <- if (identical(target, "SCZ_PGC2022")) file.path(root, "results/p11_scz_neff_audit/phbc_cache") else file.path(root, "results/p25/phbc_cache", target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target, "; found ", length(files))
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) {
    stop("incomplete or duplicated blocks for ", target)
  }
  rgarray <- array(
    NA_real_, c(length(traits), length(traits), n_block),
    dimnames = list(traits, traits, paste0("block_", seq_len(n_block)))
  )
  h2array <- matrix(
    NA_real_, n_block, length(traits),
    dimnames = list(paste0("block_", seq_len(n_block)), traits)
  )
  gcovarray <- array(
    NA_real_, c(length(traits), length(traits), n_block),
    dimnames = list(traits, traits, paste0("block_", seq_len(n_block)))
  )
  for (x in chunks) {
    if (!identical(x$target, target) || !identical(x$traits, traits)) {
      stop("cache trait-order mismatch for ", target)
    }
    idx <- x$block_ids
    rgarray[, , idx] <- x$rgarray
    h2array[idx, ] <- x$h2array
    gcovarray[, , idx] <- x$gcovarray
  }
  if (!all(is.finite(rgarray)) || !all(is.finite(h2array)) || !all(is.finite(gcovarray))) {
    stop("non-finite cache values for ", target)
  }
  list(rgarray = rgarray, h2array = h2array, gcovarray = gcovarray)
}

subset_input <- function(x, auxiliaries) {
  traits <- c(x$target, auxiliaries)
  idx <- match(traits, x$traits)
  if (anyNA(idx)) stop("missing traits for ", x$target)
  list(
    target = x$target,
    auxiliaries = auxiliaries,
    traits = traits,
    full_h2 = x$full_h2[, idx, drop = FALSE],
    full_rg = x$full_rg[idx, idx, drop = FALSE],
    h2array = x$h2array[, idx, drop = FALSE],
    rgarray = x$rgarray[idx, idx, , drop = FALSE],
    gcovarray = x$gcovarray[idx, idx, , drop = FALSE]
  )
}

uncorrected_blocks <- function(x) {
  vapply(seq_len(n_block), function(b) {
    r <- x$rgarray[1L, -1L, b]
    r_aux <- x$rgarray[-1L, -1L, b, drop = FALSE]
    dim(r_aux) <- c(length(r), length(r))
    as.numeric(crossprod(r, solve(r_aux, r)))
  }, numeric(1L))
}

fit_panel <- function(x) {
  set.seed(seed)
  result <- pleiotropyh2_cor_computing_single(
    1L, x$traits, x$full_h2, x$h2array, x$full_rg, x$rgarray, sample_rep
  )
  phbc <- as.numeric(result$percentage_h2pleio_corr)
  official_se <- as.numeric(result$percentage_h2pleio_corr_se)
  corrected_weight <- as.numeric(result$corrected_weight)
  corrected_blocks <- corrected_weight^2 * uncorrected_blocks(x)
  stable <- all(is.finite(c(phbc, official_se, corrected_weight, corrected_blocks))) &&
    corrected_weight >= 0.5 && official_se <= 0.5
  if (!stable) stop("correction stability gate failed for ", x$target)
  list(
    target = x$target,
    auxiliaries = x$auxiliaries,
    phbc = phbc,
    official_se = official_se,
    corrected_weight = corrected_weight,
    corrected_blocks = corrected_blocks
  )
}

paired_test <- function(point, block_values) {
  pseudo <- n_block * point - (n_block - 1L) * block_values
  se <- sd(pseudo) / sqrt(n_block)
  z <- point / se
  list(
    se = se,
    low = point - 1.96 * se,
    high = point + 1.96 * se,
    z = z,
    p = 2 * pnorm(-abs(z))
  )
}

base <- setNames(lapply(required_targets, function(target) {
  traits <- c(target, p25_panels[[target]])
  full_path <- if (identical(target, "SCZ_PGC2022")) file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds") else file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds")
  full <- readRDS(full_path)
  cache <- combine_chunks(target, traits)
  list(
    target = target,
    traits = traits,
    auxiliaries = p25_panels[[target]],
    full_h2 = full$h2[, traits, drop = FALSE],
    full_rg = full$rg[traits, traits, drop = FALSE],
    h2array = cache$h2array,
    rgarray = cache$rgarray,
    gcovarray = cache$gcovarray
  )
}), required_targets)

primary_records <- readRDS(file.path(root, "results/p11_scz_neff_audit/corrected_p25_final/p25_phbc_records.rds"))$full
full_records <- primary_records[required_targets]
for (target in required_targets) {
  if (!identical(full_records[[target]]$auxiliaries, p25_panels[[target]])) {
    stop("full-panel record mismatch for ", target)
  }
}

target_rows <- list()
paired_rows <- list()

append_target_row <- function(record, omitted_auxiliary, condition) {
  data.table(
    target = record$target,
    disease = ifelse(grepl("^SCZ", record$target), "SCZ", "MDD"),
    condition = condition,
    omitted_auxiliary = omitted_auxiliary,
    remaining_auxiliaries = paste(record$auxiliaries, collapse = ","),
    n_remaining = length(record$auxiliaries),
    phbc_fraction = record$phbc,
    phbc_pp = 100 * record$phbc,
    official_se_fraction = record$official_se,
    corrected_weight = record$corrected_weight,
    n_block = n_block,
    sample_rep = sample_rep
  )
}

for (target in required_targets) {
  target_rows[[length(target_rows) + 1L]] <- append_target_row(
    full_records[[target]], NA_character_, "full_panel"
  )
}

for (i in seq_len(nrow(headline_pairs))) {
  pair <- headline_pairs[i]
  left_full <- full_records[[pair$left]]
  right_full <- full_records[[pair$right]]
  if (!identical(left_full$auxiliaries, right_full$auxiliaries)) {
    stop("family panel mismatch for ", pair$comparison)
  }
  full_point <- left_full$phbc - right_full$phbc
  full_test <- paired_test(full_point, left_full$corrected_blocks - right_full$corrected_blocks)
  paired_rows[[length(paired_rows) + 1L]] <- data.table(
    comparison = pair$comparison,
    disease = pair$disease,
    condition = "full_panel",
    omitted_auxiliary = NA_character_,
    n_remaining = length(left_full$auxiliaries),
    left_phbc_pp = 100 * left_full$phbc,
    right_phbc_pp = 100 * right_full$phbc,
    delta_pp = 100 * full_point,
    paired_jackknife_se_pp = 100 * full_test$se,
    ci95_low_pp = 100 * full_test$low,
    ci95_high_pp = 100 * full_test$high,
    z = full_test$z,
    p_two_sided = full_test$p,
    deviation_from_full_delta_pp = 0,
    direction_matches_full = TRUE,
    left_corrected_weight = left_full$corrected_weight,
    right_corrected_weight = right_full$corrected_weight,
    n_block = n_block,
    sample_rep = sample_rep
  )

  for (omitted in left_full$auxiliaries) {
    remaining <- setdiff(left_full$auxiliaries, omitted)
    left <- fit_panel(subset_input(base[[pair$left]], remaining))
    right <- fit_panel(subset_input(base[[pair$right]], remaining))
    target_rows[[length(target_rows) + 1L]] <- append_target_row(left, omitted, "leave_one_auxiliary_out")
    target_rows[[length(target_rows) + 1L]] <- append_target_row(right, omitted, "leave_one_auxiliary_out")
    point <- left$phbc - right$phbc
    test <- paired_test(point, left$corrected_blocks - right$corrected_blocks)
    paired_rows[[length(paired_rows) + 1L]] <- data.table(
      comparison = pair$comparison,
      disease = pair$disease,
      condition = "leave_one_auxiliary_out",
      omitted_auxiliary = omitted,
      n_remaining = length(remaining),
      left_phbc_pp = 100 * left$phbc,
      right_phbc_pp = 100 * right$phbc,
      delta_pp = 100 * point,
      paired_jackknife_se_pp = 100 * test$se,
      ci95_low_pp = 100 * test$low,
      ci95_high_pp = 100 * test$high,
      z = test$z,
      p_two_sided = test$p,
      deviation_from_full_delta_pp = 100 * (point - full_point),
      direction_matches_full = sign(point) == sign(full_point),
      left_corrected_weight = left$corrected_weight,
      right_corrected_weight = right$corrected_weight,
      n_block = n_block,
      sample_rep = sample_rep
    )
  }
}

target_results <- unique(rbindlist(target_rows, fill = TRUE))
paired_results <- rbindlist(paired_rows, fill = TRUE)

omitted_results <- paired_results[condition == "leave_one_auxiliary_out"]
if (nrow(paired_results) != 29L || nrow(omitted_results) != 26L) {
  stop("unexpected paired result count")
}
omitted_results[, q_bh_26 := p.adjust(p_two_sided, method = "BH")]
paired_results <- merge(
  paired_results,
  omitted_results[, .(comparison, omitted_auxiliary, q_bh_26)],
  by = c("comparison", "omitted_auxiliary"), all.x = TRUE, sort = FALSE
)

summary_results <- paired_results[condition == "leave_one_auxiliary_out", .(
  n_omissions = .N,
  full_delta_pp = (delta_pp - deviation_from_full_delta_pp)[1L],
  minimum_delta_pp = min(delta_pp),
  maximum_delta_pp = max(delta_pp),
  maximum_absolute_deviation_pp = max(abs(deviation_from_full_delta_pp)),
  all_directions_match_full = all(direction_matches_full),
  n_ci_excluding_zero = sum(ci95_low_pp > 0 | ci95_high_pp < 0),
  n_bh_0_05 = sum(q_bh_26 < 0.05)
), by = .(comparison, disease)]

if (!all(summary_results$all_directions_match_full)) {
  stop("at least one leave-one-auxiliary-out estimate reversed direction")
}

setorder(target_results, disease, target, condition, omitted_auxiliary)
setorder(paired_results, disease, comparison, condition, omitted_auxiliary)
setorder(summary_results, disease, comparison)

fwrite(
  target_results,
  file.path(out_dir, "p10_leave_one_auxiliary_out_target_phbc.tsv"),
  sep = "\t"
)
fwrite(
  paired_results,
  file.path(out_dir, "p10_leave_one_auxiliary_out_paired_phbc.tsv"),
  sep = "\t"
)
fwrite(
  summary_results,
  file.path(out_dir, "p10_leave_one_auxiliary_out_summary.tsv"),
  sep = "\t"
)
saveRDS(
  list(
    parameters = list(n_block = n_block, sample_rep = sample_rep, seed = seed),
    target_results = target_results,
    paired_results = paired_results,
    summary_results = summary_results
  ),
  file.path(out_dir, "p10_leave_one_auxiliary_out_records.rds")
)

print(summary_results)
cat("P10_LEAVE_ONE_AUXILIARY_OUT_COMPLETE\n")

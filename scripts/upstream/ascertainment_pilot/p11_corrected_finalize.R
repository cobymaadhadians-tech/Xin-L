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
out_dir <- file.path(root, "results/p11_scz_neff_audit/corrected_p25_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- data.table(
  comparison = c("MDD_EHR_minus_CLIN", "MDD_FINNGEN_minus_CLIN", "MDD_QUEST_minus_CLIN", "MDD_EHR_minus_QUEST", "BD_FINNGEN_minus_CLIN", "SCZ_FINNGEN_minus_PGC"),
  disease = c("MDD", "MDD", "MDD", "MDD", "BD", "SCZ"),
  left = c("MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025", "MDD_EHR_PGC2025", "BD_FINNGEN_R13", "SCZ_FINNGEN_R13"),
  right = c("MDD_CLIN_PGC2025", "MDD_CLIN_PGC2025", "MDD_CLIN_PGC2025", "MDD_QUEST_PGC2025", "BD_CLIN_PGC4", "SCZ_PGC2022")
)

domains <- list(
  mood_psychotic = c("MDD_CLIN_PGC2025", "BD_PGC2021", "SCZ_PGC2022"),
  internalizing = c("PTSD", "ANX"),
  neurodevelopmental = c("ADHD", "ASD"),
  compulsive_eating = c("OCD_2025", "AN"),
  substance = c("AUD", "CUD_2023_EUR")
)

combine_chunks <- function(target, traits) {
  chunk_dir <- if (identical(target, "SCZ_PGC2022")) {
    file.path(root, "results/p11_scz_neff_audit/phbc_cache")
  } else file.path(root, "results/p25/phbc_cache", target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target, "; found ", length(files))
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) stop("incomplete or duplicated blocks for ", target)
  rgarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  h2array <- matrix(NA_real_, n_block, length(traits), dimnames = list(paste0("block_", seq_len(n_block)), traits))
  gcovarray <- array(NA_real_, c(length(traits), length(traits), n_block), dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  for (x in chunks) {
    if (!identical(x$target, target) || !identical(x$traits, traits)) stop("cache trait-order mismatch for ", target)
    idx <- x$block_ids
    rgarray[, , idx] <- x$rgarray
    h2array[idx, ] <- x$h2array
    gcovarray[, , idx] <- x$gcovarray
  }
  if (!all(is.finite(rgarray)) || !all(is.finite(h2array)) || !all(is.finite(gcovarray))) stop("non-finite cache values for ", target)
  list(rgarray = rgarray, h2array = h2array, gcovarray = gcovarray)
}

subset_input <- function(x, auxiliaries) {
  traits <- c(x$target, auxiliaries)
  idx <- match(traits, x$traits)
  if (anyNA(idx)) stop("missing traits in ", x$target, ": ", paste(traits[is.na(idx)], collapse = ","))
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

fit_panel <- function(x, analysis, domain = NA_character_) {
  set.seed(seed)
  result <- pleiotropyh2_cor_computing_single(1L, x$traits, x$full_h2, x$h2array, x$full_rg, x$rgarray, sample_rep)
  phbc <- as.numeric(result$percentage_h2pleio_corr)
  official_se <- as.numeric(result$percentage_h2pleio_corr_se)
  official_z <- as.numeric(result$percentage_h2pleio_corr_Z)
  weight <- as.numeric(result$corrected_weight)
  block_uncorrected <- uncorrected_blocks(x)
  corrected_blocks <- weight^2 * block_uncorrected
  pseudo <- n_block * phbc - (n_block - 1L) * corrected_blocks
  block_se <- sd(pseudo) / sqrt(n_block)
  stable <- all(is.finite(c(phbc, official_se, official_z, weight, block_se, corrected_blocks))) && weight >= 0.5 && official_se <= 0.5
  list(
    target = x$target, analysis = analysis, domain = domain, traits = x$traits,
    auxiliaries = x$auxiliaries, phbc = phbc, official_se = official_se,
    official_z = official_z, corrected_weight = weight, block_se = block_se,
    block_uncorrected = block_uncorrected, corrected_blocks = corrected_blocks,
    uncorrected_quadratic_form = uncorrected_point(x), stable = stable,
    selected_auxD = as.character(result$selected_auxD), result = result
  )
}

summary_record <- function(x) data.table(
  target = x$target,
  analysis = x$analysis,
  domain = x$domain,
  auxiliaries = paste(x$auxiliaries, collapse = ","),
  n_aux = length(x$auxiliaries),
  selected_auxD = paste(x$selected_auxD, collapse = ","),
  n_block = n_block,
  sample_rep = sample_rep,
  phbc = x$phbc,
  phbc_pp = 100 * x$phbc,
  official_se = x$official_se,
  official_z = x$official_z,
  paired_block_jackknife_se = x$block_se,
  paired_block_jackknife_se_pp = 100 * x$block_se,
  corrected_weight = x$corrected_weight,
  uncorrected_quadratic_form = x$uncorrected_quadratic_form,
  numerical_stability_gate = ifelse(x$stable, "PASS", "STOP")
)

paired_test <- function(point, blocks) {
  pseudo <- n_block * point - (n_block - 1L) * blocks
  se <- sd(pseudo) / sqrt(n_block)
  z <- point / se
  list(se = se, low = point - 1.96 * se, high = point + 1.96 * se, z = z, p = 2 * pnorm(-abs(z)))
}

# Recompute only the target whose LDSC N convention changed. Other records are
# carried forward after the auxiliary-only audit showed negligible rg changes.
recalc_targets <- "SCZ_PGC2022"
base <- setNames(lapply(recalc_targets, function(target) {
  traits <- c(target, p25_panels[[target]])
  full_path <- if (identical(target, "SCZ_PGC2022")) {
    file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds")
  } else file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds")
  full <- readRDS(full_path)
  cache <- combine_chunks(target, traits)
  list(target = target, traits = traits, auxiliaries = p25_panels[[target]],
       full_h2 = full$h2[, traits, drop = FALSE], full_rg = full$rg[traits, traits, drop = FALSE],
       h2array = cache$h2array, rgarray = cache$rgarray, gcovarray = cache$gcovarray)
}), recalc_targets)

# Full-panel primary estimates.
old_records <- readRDS(file.path(root, "results/p25/final/p25_phbc_records.rds"))
full_records <- old_records$full
full_records[recalc_targets] <- setNames(lapply(recalc_targets, function(target) fit_panel(base[[target]], "full_panel")), recalc_targets)
full_summary <- rbindlist(lapply(full_records, summary_record))
if (any(full_summary$numerical_stability_gate != "PASS")) stop("one or more full-panel PHBC estimates failed the correction stability gate")

# Six prespecified paired total-PHBC contrasts.
total_rows <- lapply(seq_len(nrow(pairs)), function(i) {
  left <- full_records[[pairs$left[[i]]]]
  right <- full_records[[pairs$right[[i]]]]
  if (!identical(left$auxiliaries, right$auxiliaries)) stop("family panel mismatch for ", pairs$comparison[[i]])
  point <- left$phbc - right$phbc
  tst <- paired_test(point, left$corrected_blocks - right$corrected_blocks)
  data.table(
    comparison = pairs$comparison[[i]], disease = pairs$disease[[i]],
    left_target = pairs$left[[i]], right_target = pairs$right[[i]],
    auxiliaries = paste(left$auxiliaries, collapse = ","), n_aux = length(left$auxiliaries),
    n_block = n_block, left_phbc = left$phbc, right_phbc = right$phbc,
    delta_fraction = point, delta_pp = 100 * point,
    paired_jackknife_se_fraction = tst$se, paired_jackknife_se_pp = 100 * tst$se,
    ci95_low_pp = 100 * tst$low, ci95_high_pp = 100 * tst$high,
    z = tst$z, p_two_sided = tst$p,
    left_corrected_weight = left$corrected_weight, right_corrected_weight = right$corrected_weight,
    method = "200-block paired delete-one-block jackknife with fixed full-data correction weights"
  )
})
total_delta <- rbindlist(total_rows)
total_delta[, `:=`(
  fdr_family = "six prespecified total-PHBC paired contrasts",
  q_bh_6 = p.adjust(p_two_sided, method = "BH")
)]
total_delta[, survives_bh_0_05 := q_bh_6 < 0.05]

# Leave-category-out estimates. Structural absence is retained as an explicit NA row.
loo_records <- old_records$leave_category_out
loo_summary_rows <- list()
for (target in names(base)) {
  for (domain in names(domains)) {
    present <- intersect(base[[target]]$auxiliaries, domains[[domain]])
    key <- paste(target, domain, sep = "__")
    if (!length(present)) {
      loo_records[[key]] <- NULL
      loo_summary_rows[[length(loo_summary_rows) + 1L]] <- data.table(
        target = target, domain = domain, removed_auxiliaries = NA_character_, n_removed = 0L,
        remaining_auxiliaries = paste(base[[target]]$auxiliaries, collapse = ","), n_remaining = length(base[[target]]$auxiliaries),
        full_phbc = full_records[[target]]$phbc, minus_category_phbc = NA_real_, leave_category_out_reduction = NA_real_,
        leave_category_out_reduction_pp = NA_real_, status = "STRUCTURAL_NA_NO_DOMAIN_MEMBER"
      )
      next
    }
    remaining <- setdiff(base[[target]]$auxiliaries, present)
    rec <- fit_panel(subset_input(base[[target]], remaining), "leave_category_out", domain)
    if (!rec$stable) stop("leave-category-out correction stability failure: ", key)
    loo_records[[key]] <- rec
    reduction <- full_records[[target]]$phbc - rec$phbc
    loo_summary_rows[[length(loo_summary_rows) + 1L]] <- data.table(
      target = target, domain = domain, removed_auxiliaries = paste(present, collapse = ","), n_removed = length(present),
      remaining_auxiliaries = paste(remaining, collapse = ","), n_remaining = length(remaining),
      full_phbc = full_records[[target]]$phbc, minus_category_phbc = rec$phbc,
      leave_category_out_reduction = reduction, leave_category_out_reduction_pp = 100 * reduction,
      status = "ESTIMABLE"
    )
  }
}
loo_summary <- rbindlist(loo_summary_rows, fill = TRUE)

# Paired difference-of-differences for each estimable family/domain combination.
loo_delta_rows <- list()
for (i in seq_len(nrow(pairs))) {
  for (domain in names(domains)) {
    left_loo <- loo_records[[paste(pairs$left[[i]], domain, sep = "__")]]
    right_loo <- loo_records[[paste(pairs$right[[i]], domain, sep = "__")]]
    if (is.null(left_loo) || is.null(right_loo)) next
    left_full <- full_records[[pairs$left[[i]]]]
    right_full <- full_records[[pairs$right[[i]]]]
    if (!identical(left_loo$auxiliaries, right_loo$auxiliaries)) stop("leave-category-out panel mismatch for ", pairs$comparison[[i]], " / ", domain)
    left_reduction <- left_full$phbc - left_loo$phbc
    right_reduction <- right_full$phbc - right_loo$phbc
    point <- left_reduction - right_reduction
    blocks <- (left_full$corrected_blocks - left_loo$corrected_blocks) - (right_full$corrected_blocks - right_loo$corrected_blocks)
    tst <- paired_test(point, blocks)
    loo_delta_rows[[length(loo_delta_rows) + 1L]] <- data.table(
      comparison = pairs$comparison[[i]], disease = pairs$disease[[i]], domain = domain,
      left_target = pairs$left[[i]], right_target = pairs$right[[i]],
      removed_auxiliaries = paste(intersect(left_full$auxiliaries, domains[[domain]]), collapse = ","),
      n_block = n_block, left_reduction_pp = 100 * left_reduction, right_reduction_pp = 100 * right_reduction,
      difference_in_reduction_pp = 100 * point, paired_jackknife_se_pp = 100 * tst$se,
      ci95_low_pp = 100 * tst$low, ci95_high_pp = 100 * tst$high,
      z = tst$z, p_two_sided = tst$p,
      method = "200-block paired difference-of-differences jackknife with panel-specific fixed correction weights"
    )
  }
}
loo_delta <- rbindlist(loo_delta_rows)
if (nrow(loo_delta) != 26L) stop("expected 26 estimable leave-category-out contrasts; found ", nrow(loo_delta))
loo_delta[, `:=`(
  fdr_family = "26 prespecified estimable leave-category-out difference-of-differences tests",
  q_bh_26 = p.adjust(p_two_sided, method = "BH")
)]
loo_delta[, survives_bh_0_05 := q_bh_26 < 0.05]

# Domain-restricted sensitivity estimates and paired contrasts.
domain_records <- old_records$domain_restricted
domain_summary_rows <- list()
for (target in names(base)) {
  for (domain in names(domains)) {
    present <- intersect(base[[target]]$auxiliaries, domains[[domain]])
    if (!length(present)) next
    rec <- fit_panel(subset_input(base[[target]], present), "domain_restricted_sensitivity", domain)
    domain_records[[paste(target, domain, sep = "__")]] <- rec
    domain_summary_rows[[length(domain_summary_rows) + 1L]] <- summary_record(rec)
  }
}
domain_summary <- rbindlist(domain_summary_rows)
domain_delta_rows <- list()
for (i in seq_len(nrow(pairs))) {
  for (domain in names(domains)) {
    left <- domain_records[[paste(pairs$left[[i]], domain, sep = "__")]]
    right <- domain_records[[paste(pairs$right[[i]], domain, sep = "__")]]
    if (is.null(left) || is.null(right)) next
    if (!identical(left$auxiliaries, right$auxiliaries)) stop("domain-restricted panel mismatch")
    point <- left$phbc - right$phbc
    tst <- paired_test(point, left$corrected_blocks - right$corrected_blocks)
    domain_delta_rows[[length(domain_delta_rows) + 1L]] <- data.table(
      comparison = pairs$comparison[[i]], disease = pairs$disease[[i]], domain = domain,
      auxiliaries = paste(left$auxiliaries, collapse = ","), n_aux = length(left$auxiliaries),
      left_phbc = left$phbc, right_phbc = right$phbc, delta_pp = 100 * point,
      paired_jackknife_se_pp = 100 * tst$se, ci95_low_pp = 100 * tst$low,
      ci95_high_pp = 100 * tst$high, z = tst$z, p_two_sided = tst$p,
      method = "domain-restricted sensitivity; 200-block paired delete-one-block jackknife"
    )
  }
}
domain_delta <- rbindlist(domain_delta_rows)
domain_delta[, q_bh_26 := p.adjust(p_two_sided, method = "BH")]
domain_delta[, survives_bh_0_05 := q_bh_26 < 0.05]

fwrite(full_summary, file.path(out_dir, "p25_full_panel_phbc.tsv"), sep = "\t")
fwrite(total_delta, file.path(out_dir, "p25_total_paired_delta_phbc.tsv"), sep = "\t")
fwrite(loo_summary, file.path(out_dir, "p25_leave_category_out_estimates.tsv"), sep = "\t")
fwrite(loo_delta, file.path(out_dir, "p25_leave_category_out_paired_did.tsv"), sep = "\t")
fwrite(domain_summary, file.path(out_dir, "p25_domain_restricted_phbc.tsv"), sep = "\t")
fwrite(domain_delta, file.path(out_dir, "p25_domain_restricted_paired_delta.tsv"), sep = "\t")
saveRDS(list(full = full_records, leave_category_out = loo_records, domain_restricted = domain_records), file.path(out_dir, "p25_phbc_records.rds"))
fwrite(data.table(
  parameter = c("n_block", "sample_rep", "seed", "target_h2_gate", "target_aux_rg_gate", "matrix_gate", "full_correction_gate", "total_FDR_family", "leave_category_out_FDR_family", "leave_category_out_interpretation", "domain_restricted_interpretation", "pleioh2g_version"),
  value = c(n_block, sample_rep, seed, "Z_h2 > 6", "rg^2 < 0.5", "finite positive-definite auxiliary rg matrix and condition number <= 1000", "corrected_weight >= 0.5 and official PHBC SE <= 0.5", "six prespecified paired contrasts", "26 prespecified estimable paired difference-of-differences tests", "category-removal estimate; categories can overlap genetically and estimates are not additive", "secondary sensitivity estimate conditional on auxiliaries within the named category", as.character(packageVersion("pleioh2g")))
), file.path(out_dir, "p25_analysis_parameters.tsv"), sep = "\t")

cat("P25_FINALIZE_COMPLETE\n")

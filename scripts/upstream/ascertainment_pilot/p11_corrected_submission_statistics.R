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
out_dir <- file.path(root, "results/p11_scz_neff_audit/corrected_p7_submission")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

intersection_aux <- c("ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR")

pairs <- data.table(
  comparison = c(
    "MDD_EHR_minus_CLIN", "MDD_FINNGEN_minus_CLIN", "MDD_QUEST_minus_CLIN",
    "MDD_EHR_minus_QUEST", "MDD_EHR_minus_FINNGEN", "MDD_FINNGEN_minus_QUEST",
    "BD_FINNGEN_minus_CLIN", "SCZ_FINNGEN_minus_PGC"
  ),
  disease = c(rep("MDD", 6L), "BD", "SCZ"),
  left = c(
    "MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025",
    "MDD_EHR_PGC2025", "MDD_EHR_PGC2025", "MDD_FINNGEN_R13",
    "BD_FINNGEN_R13", "SCZ_FINNGEN_R13"
  ),
  right = c(
    "MDD_CLIN_PGC2025", "MDD_CLIN_PGC2025", "MDD_CLIN_PGC2025",
    "MDD_QUEST_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025",
    "BD_CLIN_PGC4", "SCZ_PGC2022"
  )
)

paired_test <- function(point, blocks) {
  stopifnot(length(blocks) == n_block, all(is.finite(blocks)), is.finite(point))
  pseudo <- n_block * point - (n_block - 1L) * blocks
  se <- sd(pseudo) / sqrt(n_block)
  z <- point / se
  list(se = se, low = point - 1.96 * se, high = point + 1.96 * se,
       z = z, p = 2 * pnorm(-abs(z)))
}

pair_records <- function(records, panel_name) {
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    left <- records[[pairs$left[[i]]]]
    right <- records[[pairs$right[[i]]]]
    if (is.null(left) || is.null(right)) stop("missing record for ", pairs$comparison[[i]])
    if (!identical(left$auxiliaries, right$auxiliaries)) {
      stop("panel mismatch for ", pairs$comparison[[i]])
    }
    point <- left$phbc - right$phbc
    tst <- paired_test(point, left$corrected_blocks - right$corrected_blocks)
    data.table(
      panel = panel_name,
      comparison = pairs$comparison[[i]], disease = pairs$disease[[i]],
      left_target = pairs$left[[i]], right_target = pairs$right[[i]],
      auxiliaries = paste(left$auxiliaries, collapse = ","),
      n_aux = length(left$auxiliaries), n_block = n_block,
      left_phbc = left$phbc, right_phbc = right$phbc,
      delta_fraction = point, delta_pp = 100 * point,
      paired_jackknife_se_fraction = tst$se,
      paired_jackknife_se_pp = 100 * tst$se,
      ci95_low_pp = 100 * tst$low, ci95_high_pp = 100 * tst$high,
      z = tst$z, p_two_sided = tst$p,
      left_corrected_weight = left$corrected_weight,
      right_corrected_weight = right$corrected_weight,
      method = "200-block paired delete-one-block jackknife with fixed full-data correction weights"
    )
  })
  ans <- rbindlist(rows)
  ans[, `:=`(
    fdr_family = paste0("eight complete within-disorder psychiatric paired contrasts; ", panel_name),
    q_bh_8 = p.adjust(p_two_sided, method = "BH")
  )]
  ans[, survives_bh_0_05 := q_bh_8 < 0.05]
  ans
}

# Complete the primary disease-specific-panel family using the frozen P2.5 records.
p25_records <- readRDS(file.path(root, "results/p11_scz_neff_audit/corrected_p25_final/p25_phbc_records.rds"))
primary_delta <- pair_records(p25_records$full, "disease-specific psychiatric panel")

combine_chunks <- function(target, traits) {
  chunk_dir <- if (identical(target, "SCZ_PGC2022")) {
    file.path(root, "results/p11_scz_neff_audit/phbc_cache")
  } else file.path(root, "results/p25/phbc_cache", target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target, "; found ", length(files))
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) {
    stop("incomplete or duplicated blocks for ", target)
  }
  rgarray <- array(NA_real_, c(length(traits), length(traits), n_block),
                   dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  h2array <- matrix(NA_real_, n_block, length(traits),
                    dimnames = list(paste0("block_", seq_len(n_block)), traits))
  gcovarray <- array(NA_real_, c(length(traits), length(traits), n_block),
                     dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
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

uncorrected_blocks <- function(x) {
  vapply(seq_len(n_block), function(b) {
    r <- x$rgarray[1L, -1L, b]
    r_aux <- x$rgarray[-1L, -1L, b, drop = FALSE]
    dim(r_aux) <- c(length(r), length(r))
    as.numeric(crossprod(r, solve(r_aux, r)))
  }, numeric(1L))
}

fit_intersection <- function(target) {
  full_traits <- c(target, p25_panels[[target]])
  selected_traits <- c(target, intersection_aux)
  full_path <- if (identical(target, "SCZ_PGC2022")) {
    file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds")
  } else file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds")
  full <- readRDS(full_path)
  cache <- combine_chunks(target, full_traits)
  idx <- match(selected_traits, full_traits)
  if (anyNA(idx)) stop("intersection trait missing for ", target)
  x <- list(
    target = target,
    traits = selected_traits,
    auxiliaries = intersection_aux,
    full_h2 = full$h2[, idx, drop = FALSE],
    full_rg = full$rg[idx, idx, drop = FALSE],
    h2array = cache$h2array[, idx, drop = FALSE],
    rgarray = cache$rgarray[idx, idx, , drop = FALSE],
    gcovarray = cache$gcovarray[idx, idx, , drop = FALSE]
  )
  r_aux <- x$full_rg[-1L, -1L, drop = FALSE]
  eig <- eigen(r_aux, symmetric = TRUE, only.values = TRUE)$values
  target_aux_rg2_max <- max(x$full_rg[1L, -1L]^2)
  gate <- all(is.finite(r_aux)) && min(eig) > 0 && kappa(r_aux, exact = TRUE) <= 1000 &&
    target_aux_rg2_max < 0.5
  if (!gate) stop("intersection matrix or eligibility gate failed for ", target)

  set.seed(seed)
  result <- pleiotropyh2_cor_computing_single(
    1L, x$traits, x$full_h2, x$h2array, x$full_rg, x$rgarray, sample_rep
  )
  phbc <- as.numeric(result$percentage_h2pleio_corr)
  official_se <- as.numeric(result$percentage_h2pleio_corr_se)
  official_z <- as.numeric(result$percentage_h2pleio_corr_Z)
  weight <- as.numeric(result$corrected_weight)
  blocks <- weight^2 * uncorrected_blocks(x)
  if (!all(is.finite(c(phbc, official_se, official_z, weight, blocks))) ||
      weight < 0.5 || official_se > 0.5) {
    stop("intersection PHBC correction gate failed for ", target)
  }
  list(
    target = target, auxiliaries = intersection_aux, traits = selected_traits,
    phbc = phbc, official_se = official_se, official_z = official_z,
    corrected_weight = weight, corrected_blocks = blocks,
    min_aux_eigenvalue = min(eig), condition_number = kappa(r_aux, exact = TRUE),
    max_target_aux_rg2 = target_aux_rg2_max,
    selected_auxD = paste(as.character(result$selected_auxD), collapse = ","),
    matrix_gate = "PASS"
  )
}

intersection_records <- setNames(lapply(names(p25_panels), fit_intersection), names(p25_panels))
intersection_phbc <- rbindlist(lapply(intersection_records, function(x) data.table(
  target = x$target, auxiliaries = paste(x$auxiliaries, collapse = ","),
  n_aux = length(x$auxiliaries), n_block = n_block, sample_rep = sample_rep,
  phbc = x$phbc, phbc_pp = 100 * x$phbc,
  official_se = x$official_se, official_ci95_low = x$phbc - 1.96 * x$official_se,
  official_ci95_high = x$phbc + 1.96 * x$official_se,
  official_z = x$official_z, corrected_weight = x$corrected_weight,
  min_aux_eigenvalue = x$min_aux_eigenvalue,
  condition_number = x$condition_number,
  max_target_aux_rg2 = x$max_target_aux_rg2,
  selected_auxD = x$selected_auxD, matrix_gate = x$matrix_gate
)))
intersection_delta <- pair_records(intersection_records, "six-trait psychiatric intersection panel")

fwrite(primary_delta, file.path(out_dir, "p7_primary_eight_paired_delta.tsv"), sep = "\t")
fwrite(intersection_phbc, file.path(out_dir, "p7_intersection_phbc.tsv"), sep = "\t")
fwrite(intersection_delta, file.path(out_dir, "p7_intersection_eight_paired_delta.tsv"), sep = "\t")
fwrite(intersection_phbc[, .(
  target, auxiliaries, n_aux, min_aux_eigenvalue, condition_number,
  max_target_aux_rg2, corrected_weight, official_se, matrix_gate
)], file.path(out_dir, "p7_intersection_matrix_gate.tsv"), sep = "\t")
fwrite(data.table(
  parameter = c("n_block", "sample_rep", "seed", "intersection_auxiliaries",
                "target_h2_gate", "target_aux_rg_gate", "matrix_gate",
                "correction_gate", "primary_FDR_family", "intersection_FDR_family",
                "pleioh2g_version"),
  value = c(n_block, sample_rep, seed, paste(intersection_aux, collapse = ","),
            "Z_h2 > 6", "rg^2 < 0.5", "finite positive-definite auxiliary rg matrix and condition number <= 1000",
            "corrected_weight >= 0.5 and official PHBC SE <= 0.5",
            "eight complete within-disorder psychiatric paired contrasts",
            "eight sensitivity contrasts under the six-trait intersection panel",
            as.character(packageVersion("pleioh2g")))
), file.path(out_dir, "p7_analysis_parameters.tsv"), sep = "\t")
saveRDS(intersection_records, file.path(out_dir, "p7_intersection_records.rds"))

cat("P7_SUBMISSION_STATISTICS_COMPLETE\n")

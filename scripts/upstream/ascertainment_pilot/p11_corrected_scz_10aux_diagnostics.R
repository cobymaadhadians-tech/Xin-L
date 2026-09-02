#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages(library(data.table))
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))
n_block <- 200L
targets <- c("SCZ_PGC2022", "SCZ_FINNGEN_R13")
aux <- p25_panels[["SCZ_PGC2022"]]
stopifnot(identical(aux, p25_panels[["SCZ_FINNGEN_R13"]]), length(aux) == 10L)

combine_target <- function(target) {
  traits <- c(target, aux)
  chunk_dir <- if (identical(target, "SCZ_PGC2022")) file.path(root, "results/p11_scz_neff_audit/phbc_cache") else file.path(root, "results/p25/phbc_cache", target)
  files <- sort(list.files(chunk_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop("expected 20 chunks for ", target)
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) stop("incomplete blocks for ", target)
  rgarray <- array(NA_real_, c(length(traits), length(traits), n_block),
                   dimnames = list(traits, traits, paste0("block_", seq_len(n_block))))
  for (x in chunks) {
    if (!identical(x$traits, traits)) stop("trait mismatch for ", target)
    rgarray[, , x$block_ids] <- x$rgarray
  }
  full_path <- if (identical(target, "SCZ_PGC2022")) file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds") else file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds")
  full <- readRDS(full_path)
  list(target = target, traits = traits, auxiliaries = aux,
       full_rg = full$rg[traits, traits, drop = FALSE], rgarray = rgarray)
}

profiles <- setNames(lapply(targets, combine_target), targets)
profile_matrix <- function(x) t(vapply(seq_len(n_block), function(b) {
  x$rgarray[x$target, aux, b]
}, numeric(length(aux))))
r_pgc <- as.numeric(profiles[["SCZ_PGC2022"]]$full_rg["SCZ_PGC2022", aux])
r_fg <- as.numeric(profiles[["SCZ_FINNGEN_R13"]]$full_rg["SCZ_FINNGEN_R13", aux])
jk_pgc <- profile_matrix(profiles[["SCZ_PGC2022"]])
jk_fg <- profile_matrix(profiles[["SCZ_FINNGEN_R13"]])
colnames(jk_pgc) <- colnames(jk_fg) <- aux

delta <- r_pgc - r_fg
delta_jk <- jk_pgc - jk_fg
pseudo_delta <- sweep(delta_jk, 2L, n_block * delta,
                      function(x, y) y - (n_block - 1L) * x)
cov_delta <- cov(pseudo_delta) / n_block
rownames(cov_delta) <- colnames(cov_delta) <- aux

one <- rep(1 / length(aux), length(aux))
mean_delta <- sum(one * delta)
mean_delta_se <- sqrt(as.numeric(crossprod(one, cov_delta %*% one)))
mean_delta_z <- mean_delta / mean_delta_se
eig_delta <- eigen(cov_delta, symmetric = TRUE, only.values = TRUE)$values
if (min(eig_delta) <= 0) stop("paired profile covariance is not positive definite")
omnibus <- as.numeric(crossprod(delta, solve(cov_delta, delta)))

out_dir <- file.path(root, "results/p11_scz_neff_audit/corrected_scz_10aux_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(
  auxiliary = aux,
  SCZ_PGC_rg = r_pgc,
  SCZ_FinnGen_rg = r_fg,
  delta_PGC_minus_FinnGen = delta,
  paired_jackknife_se = sqrt(diag(cov_delta)),
  paired_z = delta / sqrt(diag(cov_delta)),
  paired_p_two_sided = 2 * pnorm(-abs(delta / sqrt(diag(cov_delta))))
), file.path(out_dir, "scz_10aux_rg_profile_paired.tsv"), sep = "\t")
write.table(cov_delta, file.path(out_dir, "scz_10aux_rg_profile_paired_cov.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
fwrite(data.table(
  n_block = n_block,
  n_aux = length(aux),
  mean_delta = mean_delta,
  mean_delta_se = mean_delta_se,
  mean_delta_z = mean_delta_z,
  mean_delta_p_two_sided = 2 * pnorm(-abs(mean_delta_z)),
  omnibus_chi2 = omnibus,
  omnibus_df = length(aux),
  omnibus_p_chisq = pchisq(omnibus, df = length(aux), lower.tail = FALSE),
  paired_cov_min_eigenvalue = min(eig_delta),
  paired_cov_max_eigenvalue = max(eig_delta),
  paired_cov_condition_number = max(eig_delta) / min(eig_delta)
), file.path(out_dir, "scz_10aux_rg_profile_tests.tsv"), sep = "\t")

# Precision-matched null using the current 10-auxiliary panel and P2.5 PHBC.
safe_cov <- function(x) {
  s <- cov(x) * nrow(x)
  e <- eigen(s, symmetric = TRUE)
  e$values[e$values < 1e-10] <- 1e-10
  e$vectors %*% (e$values * t(e$vectors))
}
rmvnorm_base <- function(n, mean, sigma) {
  e <- eigen(sigma, symmetric = TRUE)
  z <- matrix(rnorm(n * length(mean)), nrow = n)
  sweep(z %*% (e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = length(mean))), 2L, mean, "+")
}

Rdd_pgc <- profiles[["SCZ_PGC2022"]]$full_rg[aux, aux, drop = FALSE]
Rdd_fg <- profiles[["SCZ_FINNGEN_R13"]]$full_rg[aux, aux, drop = FALSE]
if (max(abs(Rdd_pgc - Rdd_fg)) > 1e-10) stop("auxiliary matrices differ between SCZ targets")
sigma_fg <- safe_cov(jk_fg)
sigma_pgc <- safe_cov(jk_pgc)
precision_center <- solve(solve(sigma_pgc) + solve(sigma_fg),
                          solve(sigma_pgc, r_pgc) + solve(sigma_fg, r_fg))
records <- readRDS(file.path(root, "results/p25/final/p25_phbc_records.rds"))$full
fg_record <- records[["SCZ_FINNGEN_R13"]]
observed_fg <- fg_record$phbc
observed_fg_weight <- fg_record$corrected_weight
inv_Rdd <- solve(Rdd_pgc)

set.seed(20260831)
n_sim <- 10000L
centers <- list(PGC = r_pgc, precision_weighted = as.numeric(precision_center))
null_rows <- lapply(names(centers), function(center_name) {
  r_sim <- rmvnorm_base(n_sim * 3L, centers[[center_name]], sigma_fg)
  valid <- apply(r_sim, 1L, function(r) {
    min(eigen(rbind(c(1, r), cbind(r, Rdd_pgc)), symmetric = TRUE, only.values = TRUE)$values) > 0
  })
  r_sim <- r_sim[valid, , drop = FALSE]
  if (nrow(r_sim) < n_sim) stop("insufficient positive-definite draws for ", center_name)
  r_sim <- r_sim[seq_len(n_sim), , drop = FALSE]
  uncorr <- rowSums((r_sim %*% inv_Rdd) * r_sim)
  corrected <- observed_fg_weight^2 * uncorr
  data.table(
    center = center_name,
    n_sim = n_sim,
    n_valid_preselection = sum(valid),
    observed_finngen_phbc = observed_fg,
    observed_finngen_weight = observed_fg_weight,
    null_mean_uncorrected = mean(uncorr),
    null_sd_uncorrected = sd(uncorr),
    null_mean_corrected_fixed_weight = mean(corrected),
    null_sd_corrected_fixed_weight = sd(corrected),
    p_precision_fixed_weight = mean(corrected <= observed_fg),
    q05_corrected_fixed_weight = unname(quantile(corrected, 0.05)),
    q50_corrected_fixed_weight = unname(quantile(corrected, 0.50)),
    q95_corrected_fixed_weight = unname(quantile(corrected, 0.95))
  )
})
fwrite(rbindlist(null_rows), file.path(out_dir, "scz_10aux_precision_matched_null.tsv"), sep = "\t")
fwrite(data.table(
  check = c("same_auxiliary_panel", "same_auxiliary_matrix", "n_block", "n_aux", "current_finngen_phbc", "current_finngen_weight"),
  value = c(TRUE, TRUE, n_block, length(aux), observed_fg, observed_fg_weight),
  detail = c(paste(aux, collapse = ","), max(abs(Rdd_pgc - Rdd_fg)), n_block,
             length(aux), observed_fg, observed_fg_weight)
), file.path(out_dir, "scz_10aux_diagnostic_parameters.tsv"), sep = "\t")
saveRDS(list(auxiliaries = aux, r_pgc = r_pgc, r_finngen = r_fg,
             delta = delta, paired_covariance = cov_delta,
             sigma_finngen = sigma_fg, sigma_pgc = sigma_pgc,
             precision_weighted_center = precision_center),
        file.path(out_dir, "scz_10aux_diagnostics.rds"))
cat("P6_SCZ_10AUX_DIAGNOSTICS_COMPLETE\n")

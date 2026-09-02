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
aux <- p25_panels[["SCZ_PGC2022"]]
if (!identical(aux, p25_panels[["SCZ_FINNGEN_R13"]])) stop("SCZ panels differ")

combine_chunks <- function(path, target) {
  files <- sort(list.files(path, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(files) != 20L) stop(target, ": expected 20 chunks")
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (!identical(sort(ids), 1:200)) stop(target, ": incomplete blocks")
  traits <- c(target, aux)
  rg <- array(NA_real_, c(11, 11, 200), dimnames = list(traits, traits, paste0("block_", 1:200)))
  h2 <- matrix(NA_real_, 200, 11, dimnames = list(paste0("block_", 1:200), traits))
  gcov <- array(NA_real_, c(11, 11, 200), dimnames = list(traits, traits, paste0("block_", 1:200)))
  for (x in chunks) {
    rg[, , x$block_ids] <- x$rgarray
    h2[x$block_ids, ] <- x$h2array
    gcov[, , x$block_ids] <- x$gcovarray
  }
  list(rgarray = rg, h2array = h2, gcovarray = gcov)
}

make_input <- function(target, corrected = FALSE) {
  if (corrected) {
    full_path <- file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds")
    cache_path <- file.path(root, "results/p11_scz_neff_audit/phbc_cache")
  } else {
    full_path <- file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds")
    cache_path <- file.path(root, "results/p25/phbc_cache", target)
  }
  full <- readRDS(full_path)
  cache <- combine_chunks(cache_path, target)
  traits <- c(target, aux)
  list(target = target, traits = traits, auxiliaries = aux,
       full_h2 = full$h2[, traits, drop = FALSE], full_rg = full$rg[traits, traits, drop = FALSE],
       h2array = cache$h2array, rgarray = cache$rgarray, gcovarray = cache$gcovarray)
}

uncorrected_blocks <- function(x) vapply(1:200, function(b) {
  r <- x$rgarray[1, -1, b]
  R <- x$rgarray[-1, -1, b, drop = FALSE]
  dim(R) <- c(10, 10)
  as.numeric(crossprod(r, solve(R, r)))
}, numeric(1))

fit_panel <- function(x) {
  set.seed(seed)
  z <- pleiotropyh2_cor_computing_single(
    1L, x$traits, x$full_h2, x$h2array, x$full_rg, x$rgarray, sample_rep
  )
  phbc <- as.numeric(z$percentage_h2pleio_corr)
  weight <- as.numeric(z$corrected_weight)
  blocks <- weight^2 * uncorrected_blocks(x)
  list(phbc = phbc, official_se = as.numeric(z$percentage_h2pleio_corr_se),
       official_z = as.numeric(z$percentage_h2pleio_corr_Z), weight = weight,
       corrected_blocks = blocks)
}

old_pgc <- fit_panel(make_input("SCZ_PGC2022", corrected = FALSE))
new_pgc <- fit_panel(make_input("SCZ_PGC2022", corrected = TRUE))
finngen <- fit_panel(make_input("SCZ_FINNGEN_R13", corrected = FALSE))

summ <- rbindlist(list(
  data.table(source = "PGC_original_NEFF", phbc = old_pgc$phbc, official_se = old_pgc$official_se,
             official_z = old_pgc$official_z, corrected_weight = old_pgc$weight),
  data.table(source = "PGC_corrected_2xNEFF", phbc = new_pgc$phbc, official_se = new_pgc$official_se,
             official_z = new_pgc$official_z, corrected_weight = new_pgc$weight),
  data.table(source = "FinnGen", phbc = finngen$phbc, official_se = finngen$official_se,
             official_z = finngen$official_z, corrected_weight = finngen$weight)
))

paired <- function(left, right, label) {
  point <- left$phbc - right$phbc
  block_delta <- left$corrected_blocks - right$corrected_blocks
  pseudo <- n_block * point - (n_block - 1) * block_delta
  se <- sd(pseudo) / sqrt(n_block)
  data.table(comparison = label, delta_fraction = point, delta_pp = 100 * point,
             paired_se_fraction = se, paired_se_pp = 100 * se,
             ci95_low_pp = 100 * (point - 1.96 * se), ci95_high_pp = 100 * (point + 1.96 * se),
             z = point / se, p_two_sided = 2 * pnorm(-abs(point / se)))
}

contrasts <- rbindlist(list(
  paired(finngen, old_pgc, "FinnGen_minus_PGC_original_NEFF"),
  paired(finngen, new_pgc, "FinnGen_minus_PGC_corrected_2xNEFF")
))

old_full <- readRDS(file.path(root, "results/p25/full_matrix/SCZ_PGC2022/ldsc_alltraits_observed.rds"))
new_full <- readRDS(file.path(root, "results/p11_scz_neff_audit/full_matrix/ldsc_alltraits_observed_Nx2.rds"))
matrix_cmp <- data.table(
  metric = c("max_abs_rg_change", "SCZ_h2_ratio_new_old", "median_SCZ_gcov_ratio_new_old"),
  value = c(max(abs(new_full$rg - old_full$rg)),
            new_full$h2[1, "SCZ_PGC2022"] / old_full$h2[1, "SCZ_PGC2022"],
            median(new_full$gcov["SCZ_PGC2022", ] / old_full$gcov["SCZ_PGC2022", ]))
)

out <- file.path(root, "results/p11_scz_neff_audit")
fwrite(summ, file.path(out, "phbc_old_vs_corrected.tsv"), sep = "\t")
fwrite(contrasts, file.path(out, "paired_contrast_old_vs_corrected.tsv"), sep = "\t")
fwrite(matrix_cmp, file.path(out, "matrix_scale_old_vs_corrected.tsv"), sep = "\t")
writeLines(c("comparison_complete=TRUE", paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))),
           file.path(out, "p11_neff_status.txt"))

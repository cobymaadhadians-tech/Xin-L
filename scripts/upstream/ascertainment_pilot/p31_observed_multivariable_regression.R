#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 12)

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
outdir <- file.path(root, "results", "p31_scz_observed_multivariable_regression")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source_ids <- c(PGC = "SCZ_PGC2022", FinnGen = "SCZ_FINNGEN_R13")
primary_predictors <- c("BD_PGC2021", "MDD_CLIN_PGC2025", "ADHD", "OCD_2025", "PTSD", "AUD")
sensitivity_predictors <- c(
  "MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD",
  "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"
)
B <- 200L

load_p25_source <- function(source_id) {
  full_path <- file.path(root, "results", "p25", "full_matrix", source_id, "ldsc_alltraits_observed.rds")
  full <- readRDS(full_path)
  chunk_paths <- file.path(
    root, "results", "p25", "phbc_cache", source_id,
    sprintf("chunk_%03d.rds", 0:19)
  )
  chunks <- lapply(chunk_paths, readRDS)
  block_ids <- unlist(lapply(chunks, function(x) x$block_ids), use.names = FALSE)
  if (length(block_ids) != B || !identical(as.integer(block_ids), seq_len(B))) {
    stop("Expected ordered block IDs 1:200 for ", source_id)
  }
  rgarray <- abind::abind(lapply(chunks, function(x) x$rgarray), along = 3)
  gcovarray <- abind::abind(lapply(chunks, function(x) x$gcovarray), along = 3)
  h2array <- do.call(rbind, lapply(chunks, function(x) x$h2array))
  if (dim(rgarray)[3] != B || dim(gcovarray)[3] != B || nrow(h2array) != B) {
    stop("Incomplete block arrays for ", source_id)
  }
  list(full = full, rgarray = rgarray, gcovarray = gcovarray, h2array = h2array)
}

load_corrected_pgc <- function() {
  base <- file.path(root, "results", "p11_scz_neff_audit")
  full <- readRDS(file.path(base, "full_matrix", "ldsc_alltraits_observed_Nx2.rds"))
  chunk_paths <- file.path(base, "phbc_cache", sprintf("chunk_%03d.rds", 0:19))
  chunks <- lapply(chunk_paths, readRDS)
  block_ids <- unlist(lapply(chunks, function(x) x$block_ids), use.names = FALSE)
  if (length(block_ids) != B || !identical(as.integer(block_ids), seq_len(B))) {
    stop("Expected ordered corrected PGC block IDs 1:200")
  }
  rgarray <- abind::abind(lapply(chunks, function(x) x$rgarray), along = 3)
  gcovarray <- abind::abind(lapply(chunks, function(x) x$gcovarray), along = 3)
  h2array <- do.call(rbind, lapply(chunks, function(x) x$h2array))
  list(full = full, rgarray = rgarray, gcovarray = gcovarray, h2array = h2array)
}

matrix_metrics <- function(Rx) {
  ev <- eigen(Rx, symmetric = TRUE, only.values = TRUE)$values
  data.frame(
    min_eigenvalue = min(ev),
    max_eigenvalue = max(ev),
    positive_definite = min(ev) > 0,
    condition_number = kappa(Rx, exact = TRUE),
    max_vif = max(diag(solve(Rx)))
  )
}

fit_standardized <- function(R, target, predictors) {
  idx <- match(c(target, predictors), colnames(R))
  if (anyNA(idx)) stop("Missing trait(s): ", paste(c(target, predictors)[is.na(idx)], collapse = ", "))
  Rt <- R[idx, idx, drop = FALSE]
  Rx <- Rt[-1, -1, drop = FALSE]
  rxy <- Rt[-1, 1, drop = FALSE]
  beta <- as.numeric(solve(Rx, rxy))
  names(beta) <- predictors
  R2 <- as.numeric(t(rxy) %*% solve(Rx, rxy))
  list(beta = beta, R2 = R2, metrics = matrix_metrics(Rx))
}

jk_inference <- function(full_estimate, delete_one) {
  full_estimate <- unname(full_estimate)
  delete_one <- unname(delete_one)
  n <- length(delete_one)
  center <- mean(delete_one)
  se <- sqrt((n - 1) / n * sum((delete_one - center)^2))
  z <- full_estimate / se
  p <- 2 * pnorm(-abs(z))
  c(
    estimate = full_estimate,
    jackknife_mean = center,
    se = se,
    ci_low = full_estimate - qnorm(0.975) * se,
    ci_high = full_estimate + qnorm(0.975) * se,
    z = z,
    p = p,
    min_delete_one = min(delete_one),
    max_delete_one = max(delete_one),
    max_abs_deletion_deviation = max(abs(delete_one - full_estimate)),
    max_influence_block = which.max(abs(delete_one - full_estimate)),
    direction_reversals = sum(sign(delete_one) != sign(full_estimate))
  )
}

sources <- list(
  PGC = load_corrected_pgc(),
  FinnGen = load_p25_source(source_ids[["FinnGen"]])
)

target_name <- function(x) colnames(x$full$rg)[1]
targets <- vapply(sources, target_name, character(1))

# The auxiliary-only submatrices should be identical because all predictors are
# the same GWASs and only the target SCZ GWAS changes.
aux_names <- sensitivity_predictors
full_aux_diff <- max(abs(
  sources$PGC$full$rg[aux_names, aux_names] -
    sources$FinnGen$full$rg[aux_names, aux_names]
))
block_aux_diff <- max(vapply(seq_len(B), function(b) {
  max(abs(
    sources$PGC$rgarray[aux_names, aux_names, b] -
      sources$FinnGen$rgarray[aux_names, aux_names, b]
  ))
}, numeric(1)))
if (full_aux_diff > 1e-10 || block_aux_diff > 1e-10) {
  stop("Predictor correlation matrices differ across target-source runs")
}

run_model <- function(model_name, predictors) {
  full_fits <- lapply(names(sources), function(s) {
    fit_standardized(sources[[s]]$full$rg, targets[[s]], predictors)
  })
  names(full_fits) <- names(sources)

  block_fits <- lapply(names(sources), function(s) {
    beta <- matrix(NA_real_, nrow = B, ncol = length(predictors),
                   dimnames = list(NULL, predictors))
    R2 <- numeric(B)
    metrics <- vector("list", B)
    for (b in seq_len(B)) {
      f <- fit_standardized(sources[[s]]$rgarray[, , b], targets[[s]], predictors)
      beta[b, ] <- f$beta
      R2[b] <- f$R2
      metrics[[b]] <- f$metrics
    }
    list(beta = beta, R2 = R2, metrics = do.call(rbind, metrics))
  })
  names(block_fits) <- names(sources)

  delta_beta_full <- full_fits$FinnGen$beta - full_fits$PGC$beta
  delta_beta_block <- block_fits$FinnGen$beta - block_fits$PGC$beta
  coef_rows <- do.call(rbind, lapply(seq_along(predictors), function(j) {
    inf <- jk_inference(delta_beta_full[j], delta_beta_block[, j])
    data.frame(
      model = model_name,
      predictor = predictors[j],
      beta_PGC = full_fits$PGC$beta[j],
      beta_FinnGen = full_fits$FinnGen$beta[j],
      delta_beta_FinnGen_minus_PGC = inf[["estimate"]],
      paired_jackknife_SE = inf[["se"]],
      CI95_low = inf[["ci_low"]],
      CI95_high = inf[["ci_high"]],
      z = inf[["z"]],
      P = inf[["p"]],
      min_delete_one_delta = inf[["min_delete_one"]],
      max_delete_one_delta = inf[["max_delete_one"]],
      max_abs_deletion_deviation = inf[["max_abs_deletion_deviation"]],
      max_influence_block = as.integer(inf[["max_influence_block"]]),
      direction_reversals_200 = as.integer(inf[["direction_reversals"]])
    )
  }))
  coef_rows$FDR_BH_within_model <- p.adjust(coef_rows$P, method = "BH")

  delta_R2_full <- full_fits$FinnGen$R2 - full_fits$PGC$R2
  delta_R2_block <- block_fits$FinnGen$R2 - block_fits$PGC$R2
  r2_inf <- jk_inference(delta_R2_full, delta_R2_block)
  r2_row <- data.frame(
    model = model_name,
    R2_PGC = full_fits$PGC$R2,
    R2_FinnGen = full_fits$FinnGen$R2,
    delta_R2_FinnGen_minus_PGC = r2_inf[["estimate"]],
    paired_jackknife_SE = r2_inf[["se"]],
    CI95_low = r2_inf[["ci_low"]],
    CI95_high = r2_inf[["ci_high"]],
    z = r2_inf[["z"]],
    P = r2_inf[["p"]],
    min_delete_one_delta = r2_inf[["min_delete_one"]],
    max_delete_one_delta = r2_inf[["max_delete_one"]],
    max_abs_deletion_deviation = r2_inf[["max_abs_deletion_deviation"]],
    max_influence_block = as.integer(r2_inf[["max_influence_block"]]),
    direction_reversals_200 = as.integer(r2_inf[["direction_reversals"]])
  )

  diag_rows <- do.call(rbind, lapply(names(sources), function(s) {
    full_m <- full_fits[[s]]$metrics
    bm <- block_fits[[s]]$metrics
    data.frame(
      model = model_name,
      source = s,
      full_min_eigenvalue = full_m$min_eigenvalue,
      full_positive_definite = full_m$positive_definite,
      full_condition_number = full_m$condition_number,
      full_max_VIF = full_m$max_vif,
      min_block_eigenvalue = min(bm$min_eigenvalue),
      max_block_condition_number = max(bm$condition_number),
      max_block_VIF = max(bm$max_vif),
      non_positive_definite_blocks = sum(!bm$positive_definite),
      inadmissible_R2_blocks = sum(block_fits[[s]]$R2 < 0 | block_fits[[s]]$R2 > 1)
    )
  }))

  block_table <- data.frame(
    model = model_name,
    block = seq_len(B),
    R2_PGC = block_fits$PGC$R2,
    R2_FinnGen = block_fits$FinnGen$R2,
    delta_R2 = delta_R2_block
  )
  for (j in seq_along(predictors)) {
    nm <- predictors[j]
    block_table[[paste0("beta_PGC__", nm)]] <- block_fits$PGC$beta[, j]
    block_table[[paste0("beta_FinnGen__", nm)]] <- block_fits$FinnGen$beta[, j]
    block_table[[paste0("delta_beta__", nm)]] <- delta_beta_block[, j]
  }

  list(coefficients = coef_rows, R2 = r2_row, diagnostics = diag_rows,
       blocks = block_table, full = full_fits)
}

primary <- run_model("primary_six", primary_predictors)
sensitivity <- run_model("sensitivity_all_ten", sensitivity_predictors)

coef_all <- rbind(primary$coefficients, sensitivity$coefficients)
r2_all <- rbind(primary$R2, sensitivity$R2)
diag_all <- rbind(primary$diagnostics, sensitivity$diagnostics)
validation <- data.frame(
  check = c("full_auxiliary_matrix_max_abs_difference", "block_auxiliary_matrix_max_abs_difference"),
  value = c(full_aux_diff, block_aux_diff),
  passed = c(full_aux_diff <= 1e-10, block_aux_diff <= 1e-10)
)

write.table(coef_all, file.path(outdir, "conditional_coefficients.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(r2_all, file.path(outdir, "joint_genetic_R2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(diag_all, file.path(outdir, "predictor_matrix_diagnostics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(validation, file.path(outdir, "matched_predictor_validation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(primary$blocks, file.path(outdir, "primary_six_delete_one_blocks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sensitivity$blocks, file.path(outdir, "sensitivity_ten_delete_one_blocks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(list(primary = primary, sensitivity = sensitivity, validation = validation),
        file.path(outdir, "observed_multivariable_regression_complete.rds"))

cat("Completed observed-variable multivariable genetic regression\n")
print(coef_all)
print(r2_all)
print(diag_all)
print(validation)

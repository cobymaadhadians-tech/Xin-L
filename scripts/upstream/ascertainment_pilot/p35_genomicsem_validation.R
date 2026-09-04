#!/usr/bin/env Rscript

# D5: Published-model Genomic SEM transport validation.
#
# The published four-factor topology is represented in Figure 1b of
# Grotzinger et al. (2022). The hierarchical fallback is represented in
# Figure 1c. Supplementary Table 52 is an alcohol-use GWAS catalogue table;
# it is not the source of the factor topology.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

root <- normalizePath(
  Sys.getenv("ANALYSIS_ROOT", unset = getwd()),
  mustWork = FALSE
)

target <- "SCZ_PGC2022"
source_dir <- file.path(
  root, "results", "p30_scz_downstream_validation", target
)
out_dir <- file.path(
  root, "results", "p35_genomicsem_validation", target
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(source_dir, "model_covstruc.rds"))) {
  stop("Missing independently estimated p30 S/V object: ", source_dir)
}

# Figure 1b: correlated factors model, including the published cross-loadings.
four_factor_model <- paste(
  "F_COMP =~ AN + OCD + TS",
  "F_PSY =~ SCZ + BIP + ALCH",
  "F_ND =~ TS + ALCH + ADHD + ASD + PTSD + MDD",
  "F_INT =~ ALCH + PTSD + MDD + ANX",
  "F_COMP ~~ F_PSY + F_ND + F_INT",
  "F_PSY ~~ F_ND + F_INT",
  "F_ND ~~ F_INT",
  sep = "\n"
)

# Figure 1c: the same first-order indicator structure with a higher-order
# p factor and no residual covariances among the first-order factors.
hierarchical_model <- paste(
  "F_COMP =~ AN + OCD + TS",
  "F_PSY =~ SCZ + BIP + ALCH",
  "F_ND =~ TS + ALCH + ADHD + ASD + PTSD + MDD",
  "F_INT =~ ALCH + PTSD + MDD + ANX",
  "P =~ F_COMP + F_PSY + F_ND + F_INT",
  sep = "\n"
)

covstruc_object <- readRDS(file.path(source_dir, "model_covstruc.rds"))
covstruc <- covstruc_object$covstruc
if (is.null(covstruc) || length(covstruc) != 2L) {
  stop("The p30 object does not contain a two-element covstruc object")
}

observed <- c("SCZ", "BIP", "MDD", "ANX", "PTSD", "ADHD", "ASD", "AN", "OCD", "TS", "ALCH")
if (!identical(sort(rownames(covstruc$S)), sort(observed))) {
  stop("Unexpected observed-trait universe in covstruc")
}

run_model <- function(model, model_id, latent_names) {
  captured_warnings <- character()
  fit <- withCallingHandlers(
    usermodel(
      covstruc = covstruc,
      estimation = "DWLS",
      model = model,
      CFIcalc = TRUE,
      std.lv = TRUE,
      imp_cov = TRUE,
      fix_resid = FALSE
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (is.null(fit) || is.null(fit$results) || !nrow(as.data.frame(fit$results))) {
    stop(model_id, ": GenomicSEM returned no interpretable parameter table")
  }

  params <- as.data.table(fit$results)
  fit_metrics <- as.data.table(fit$modelfit)
  std_col <- if ("STD_All" %in% names(params)) "STD_All" else "STD_Genotype"
  if (!std_col %in% names(params)) stop(model_id, ": no standardized estimate column")

  residuals <- params[op == "~~" & lhs == rhs & lhs %chin% observed]
  loadings <- params[op == "=~"]
  latent_variances <- params[
    op == "~~" & lhs == rhs & lhs %chin% latent_names
  ]
  factor_correlations <- params[
    op == "~~" & lhs != rhs & lhs %chin% latent_names & rhs %chin% latent_names
  ]

  cfi <- as.numeric(fit_metrics$CFI[[1L]])
  srmr <- as.numeric(fit_metrics$SRMR[[1L]])
  checks <- data.table(
    target_source = target,
    model_id = model_id,
    returned_parameter_table = TRUE,
    converged = TRUE,
    CFI = cfi,
    SRMR = srmr,
    cfi_pass = is.finite(cfi) && cfi >= 0.95,
    srmr_pass = is.finite(srmr) && srmr <= 0.08,
    residuals_nonnegative = nrow(residuals) == length(observed) && all(residuals[[std_col]] >= 0),
    latent_variances_positive = nrow(latent_variances) == length(latent_names) && all(latent_variances[[std_col]] > 0),
    factor_correlations_in_bounds = if (length(latent_names) == 5L) {
      nrow(factor_correlations) == 0L
    } else {
      nrow(factor_correlations) == 6L && all(abs(factor_correlations[[std_col]]) < 1)
    },
    standardized_loadings_in_bounds = all(is.finite(loadings[[std_col]]) & abs(loadings[[std_col]]) <= 1),
    published_loading_directions_consistent = all(loadings[[std_col]] > 0),
    warning_count = length(unique(captured_warnings))
  )
  checks[, no_heywood := residuals_nonnegative & latent_variances_positive &
    factor_correlations_in_bounds & standardized_loadings_in_bounds]
  checks[, gate_pass := converged & cfi_pass & srmr_pass & no_heywood &
    published_loading_directions_consistent & warning_count == 0L]

  params[, `:=`(target_source = target, model_id = model_id, standardized_column = std_col)]
  fit_metrics[, `:=`(target_source = target, model_id = model_id)]
  residuals[, `:=`(target_source = target, model_id = model_id, standardized_column = std_col)]
  loadings[, `:=`(target_source = target, model_id = model_id, standardized_column = std_col)]
  latent_variances[, `:=`(target_source = target, model_id = model_id, standardized_column = std_col)]
  factor_correlations[, `:=`(target_source = target, model_id = model_id, standardized_column = std_col)]

  saveRDS(fit, file.path(out_dir, paste0(model_id, "_fit.rds")))
  fwrite(params, file.path(out_dir, paste0(model_id, "_parameters.tsv")), sep = "\t")
  fwrite(fit_metrics, file.path(out_dir, paste0(model_id, "_fit_metrics.tsv")), sep = "\t")
  fwrite(residuals, file.path(out_dir, paste0(model_id, "_indicator_residuals.tsv")), sep = "\t")
  fwrite(loadings, file.path(out_dir, paste0(model_id, "_standardized_loadings.tsv")), sep = "\t")
  fwrite(latent_variances, file.path(out_dir, paste0(model_id, "_latent_variances.tsv")), sep = "\t")
  fwrite(factor_correlations, file.path(out_dir, paste0(model_id, "_factor_correlations.tsv")), sep = "\t")
  writeLines(unique(captured_warnings), file.path(out_dir, paste0(model_id, "_warnings.txt")))
  list(fit = fit, checks = checks)
}

writeLines(four_factor_model, file.path(out_dir, "published_four_factor_model.txt"))
writeLines(hierarchical_model, file.path(out_dir, "published_hierarchical_model.txt"))

# The four-factor PGC gate was already run in p30 and failed. This script runs
# only the single prespecified hierarchical fallback; FinnGen is not started.
hier <- run_model(
  hierarchical_model,
  "hierarchical_fallback",
  c("F_COMP", "F_PSY", "F_ND", "F_INT", "P")
)

fwrite(hier$checks, file.path(out_dir, "hierarchical_fallback_gate.tsv"), sep = "\t")
writeLines(c(
  "analysis = D5 published-model Genomic SEM transport validation",
  "source_model = Grotzinger et al. 2022 Nature Genetics Figure 1b/1c",
  "supplementary_table_52_correction = Table 52 is an alcohol-use GWAS catalogue table, not the factor topology",
  "four_factor_pgc_gate = FAIL in existing p30 validation; hierarchical fallback run",
  "fingen_run = NOT_STARTED because PGC reference gate did not pass",
  paste0("generated = ", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(out_dir, "D5_STATUS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat(
  "D5_HIERARCHICAL_FALLBACK\t", target, "\t",
  if (hier$checks$gate_pass) "PASS" else "FAIL", "\n", sep = ""
)

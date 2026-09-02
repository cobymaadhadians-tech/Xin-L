#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out <- file.path(root, "results", "p11_genomicsem_scz_swap")

pfiles <- sort(list.files(out, pattern = "^block_parameters_[0-9]{2}\\.tsv$", full.names = TRUE))
ffiles <- sort(list.files(out, pattern = "^block_fit_[0-9]{2}\\.tsv$", full.names = TRUE))
if (length(pfiles) != 20L || length(ffiles) != 20L) stop("Expected 20 parameter and 20 fit files")

bp <- rbindlist(lapply(pfiles, fread), fill = TRUE)
bf <- rbindlist(lapply(ffiles, fread), fill = TRUE)
if (bp[converged == TRUE, uniqueN(block)] != 200L) stop("Not all blocks converged")

fullp <- fread(file.path(out, "full_model_parameters.tsv"))
fullf <- fread(file.path(out, "full_model_fit.tsv"))

key_params <- data.table(
  lhs = c("F_PSY", "SCZ", "F_PSY"),
  op = c("=~", "~~", "~~"),
  rhs = c("SCZ", "SCZ", "F_SUB"),
  quantity = c("SCZ_standardized_loading", "SCZ_standardized_residual_variance", "factor_correlation")
)

extract_pairs <- function(dat, value, model_filter = NULL) {
  x <- copy(dat)
  if (!is.null(model_filter) && "model_id" %in% names(x)) x <- x[model_id == model_filter]
  x <- merge(x, key_params, by = c("lhs", "op", "rhs"))
  dcast(x, block + quantity ~ target_source, value.var = value)
}

blocks <- extract_pairs(bp[converged == TRUE], "STD_All")
blocks[, delta_finngen_minus_pgc := SCZ_FINNGEN_R13 - SCZ_PGC2022]

full_primary <- merge(fullp[model_id == "primary_two_factor"], key_params, by = c("lhs", "op", "rhs"))
fullwide <- dcast(full_primary, quantity ~ target_source, value.var = "STD_All")
fullwide[, delta_finngen_minus_pgc := SCZ_FINNGEN_R13 - SCZ_PGC2022]

summary <- blocks[, {
  B <- .N
  m <- mean(delta_finngen_minus_pgc)
  se <- sqrt((B - 1) / B * sum((delta_finngen_minus_pgc - m)^2))
  .(n_blocks = B, jackknife_mean_delta = m, paired_jackknife_se = se)
}, by = quantity]
summary <- merge(fullwide, summary, by = "quantity")
summary[, `:=`(
  ci95_low = delta_finngen_minus_pgc - 1.96 * paired_jackknife_se,
  ci95_high = delta_finngen_minus_pgc + 1.96 * paired_jackknife_se,
  z = delta_finngen_minus_pgc / paired_jackknife_se,
  p_two_sided = 2 * pnorm(-abs(delta_finngen_minus_pgc / paired_jackknife_se))
)]
fwrite(summary, file.path(out, "primary_parameter_swap_comparison.tsv"), sep = "\t")

fit_metrics <- intersect(c("chisq", "df", "p_chisq", "AIC", "CFI", "SRMR"), names(bf))
fit_long <- melt(bf, id.vars = c("target_source", "block"), measure.vars = fit_metrics,
                 variable.name = "metric", value.name = "value")
fit_wide <- dcast(fit_long, block + metric ~ target_source, value.var = "value")
fit_wide[, delta_finngen_minus_pgc := SCZ_FINNGEN_R13 - SCZ_PGC2022]
fit_jk <- fit_wide[, {
  B <- .N
  m <- mean(delta_finngen_minus_pgc)
  se <- sqrt((B - 1) / B * sum((delta_finngen_minus_pgc - m)^2))
  .(n_blocks = B, jackknife_mean_delta = m, paired_jackknife_se = se)
}, by = metric]

full_fit_long <- melt(fullf[model_id == "primary_two_factor"],
                      id.vars = c("target_source", "model_id"), measure.vars = fit_metrics,
                      variable.name = "metric", value.name = "value")
full_fit_wide <- dcast(full_fit_long, metric ~ target_source, value.var = "value")
full_fit_wide[, delta_finngen_minus_pgc := SCZ_FINNGEN_R13 - SCZ_PGC2022]
fit_summary <- merge(full_fit_wide, fit_jk, by = "metric", all.x = TRUE)
fwrite(fit_summary, file.path(out, "primary_fit_swap_comparison.tsv"), sep = "\t")

general_params <- fullp[model_id == "general_one_factor" &
                          ((lhs == "F_P" & op == "=~" & rhs == "SCZ") |
                           (lhs == "SCZ" & op == "~~" & rhs == "SCZ"))]
general_params[, quantity := fifelse(op == "=~", "SCZ_standardized_loading", "SCZ_standardized_residual_variance")]
general_wide <- dcast(general_params, quantity ~ target_source, value.var = "STD_All")
general_wide[, delta_finngen_minus_pgc := SCZ_FINNGEN_R13 - SCZ_PGC2022]
fwrite(general_wide, file.path(out, "general_factor_descriptive_swap.tsv"), sep = "\t")

writeLines(c(
  "analysis_complete=TRUE",
  paste0("primary_blocks=", uniqueN(blocks$block)),
  paste0("all_primary_block_fits_converged=", all(bp$converged %in% TRUE)),
  paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(out, "p11_status.txt"))

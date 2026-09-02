#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages(library(data.table))
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

rows <- vector("list", length(p25_panels))
for (i in seq_along(p25_panels)) {
  target <- names(p25_panels)[[i]]
  aux <- p25_panels[[i]]
  out_dir <- file.path(root, "results/p25/full_matrix", target)
  rds <- file.path(out_dir, "ldsc_alltraits_observed.rds")
  if (!file.exists(rds)) stop("missing full matrix RDS for ", target)
  x <- readRDS(rds)
  write.table(x$rg, file.path(out_dir, "ldsc_rg_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(x$h2Z, file.path(out_dir, "ldsc_h2z.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  aux_rg <- x$rg[aux, aux, drop = FALSE]
  finite_matrix <- all(is.finite(aux_rg))
  eigenvalues <- if (finite_matrix) eigen(aux_rg, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  min_eigenvalue <- if (finite_matrix) min(eigenvalues) else NA_real_
  condition_number <- if (finite_matrix && min_eigenvalue > 0) max(eigenvalues) / min_eigenvalue else Inf
  matrix_gate <- if (!finite_matrix || min_eigenvalue <= 0 || condition_number > 1e4) "STOP" else if (condition_number > 1e3) "CONDITIONAL" else "PASS"
  rows[[i]] <- data.table(
    target = target,
    auxiliaries = paste(aux, collapse = ","),
    n_aux = length(aux),
    finite_matrix = finite_matrix,
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = if (finite_matrix) max(eigenvalues) else NA_real_,
    condition_number = condition_number,
    matrix_gate = matrix_gate
  )
  fwrite(rows[[i]], file.path(out_dir, "matrix_gate.tsv"), sep = "\t")
}
summary <- rbindlist(rows)
fwrite(summary, file.path(root, "results/p25/full_matrix/matrix_gate_summary.tsv"), sep = "\t")
if (any(summary$matrix_gate == "STOP")) stop("at least one full matrix failed stability gate")
cat("P25_MATRIX_POSTPROCESS_COMPLETE\n")

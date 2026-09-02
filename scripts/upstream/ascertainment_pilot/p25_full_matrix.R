#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("target key required")
target_key <- args[[1L]]
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))
if (!target_key %in% names(p25_panels)) stop("unsupported target: ", target_key)

aux <- p25_panels[[target_key]]
traits <- c(target_key, aux)
all_paths <- p25_sumstats_paths(root)
paths <- all_paths[traits]
if (anyNA(paths) || !all(file.exists(paths))) stop("missing input for ", target_key)
munged <- lapply(paths, function(path) fread(path, showProgress = FALSE))
names(munged) <- traits

ld_path <- Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr"))
wld_path <- Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC"))
all_rg <- Cal_rg_h2g_alltraits(
  phenotype = traits,
  munged_sumstats = munged,
  ld_path = ld_path,
  wld_path = wld_path,
  sample_prev = NULL,
  population_prev = NULL
)

aux_rg <- all_rg$rg[aux, aux, drop = FALSE]
finite_matrix <- all(is.finite(aux_rg))
eigenvalues <- if (finite_matrix) eigen(aux_rg, symmetric = TRUE, only.values = TRUE)$values else NA_real_
min_eigenvalue <- if (finite_matrix) min(eigenvalues) else NA_real_
condition_number <- if (finite_matrix && min_eigenvalue > 0) max(eigenvalues) / min_eigenvalue else Inf
matrix_gate <- if (!finite_matrix || min_eigenvalue <= 0 || condition_number > 1e4) "STOP" else if (condition_number > 1e3) "CONDITIONAL" else "PASS"

out_dir <- file.path(root, "results/p25/full_matrix", target_key)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(all_rg, file.path(out_dir, "ldsc_alltraits_observed.rds"))
write.table(all_rg$rg, file.path(out_dir, "ldsc_rg_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(all_rg$h2Z, file.path(out_dir, "ldsc_h2z.tsv"), sep = "\t", quote = FALSE, col.names = NA)
fwrite(data.table(
  target = target_key,
  auxiliaries = paste(aux, collapse = ","),
  n_aux = length(aux),
  finite_matrix = finite_matrix,
  min_eigenvalue = min_eigenvalue,
  max_eigenvalue = if (finite_matrix) max(eigenvalues) else NA_real_,
  condition_number = condition_number,
  matrix_gate = matrix_gate
), file.path(out_dir, "matrix_gate.tsv"), sep = "\t")
if (matrix_gate == "STOP") stop("matrix stability STOP for ", target_key)
cat("P25_FULL_MATRIX_COMPLETE\t", target_key, "\t", matrix_gate, "\n", sep = "")

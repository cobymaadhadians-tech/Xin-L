#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
base_dir <- file.path(root, "results/p5/neurological_auxiliary")
rg_dir <- file.path(base_dir, "ldsc_rg")
out_dir <- file.path(base_dir, "matrix_gate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

targets <- c(
  "AD_GCST90704646_MAIN", "AD_GCST90704647_NOPROXY", "AD_GCST90704648_NOBIOBANK",
  "EPILEPSY_ILAE2023_EUR", "EPILEPSY_EHR_META_2OF3", "EPILEPSY_EHR_META_3OF3"
)
aux <- c(
  "NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641",
  "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387"
)

read_rg <- function(left, right) {
  path <- file.path(rg_dir, paste0(left, "__", right, ".log"))
  if (!file.exists(path)) stop("missing rg log: ", path)
  x <- readLines(path, warn = FALSE)
  line <- grep("^Genetic Correlation:", x, value = TRUE)
  if (length(line) != 1L) stop("rg estimate unavailable: ", path)
  m <- regexec("Genetic Correlation: ([^ ]+) \\(([^)]+)\\)", line)
  z <- regmatches(line, m)[[1L]]
  if (length(z) != 3L) stop("cannot parse rg: ", path)
  c(rg = as.numeric(z[[2L]]), se = as.numeric(z[[3L]]))
}

aux_rg <- diag(length(aux))
dimnames(aux_rg) <- list(aux, aux)
aux_se <- matrix(0, length(aux), length(aux), dimnames = list(aux, aux))
pair_rows <- list()
for (i in seq_len(length(aux) - 1L)) {
  for (j in (i + 1L):length(aux)) {
    z <- read_rg(aux[[i]], aux[[j]])
    aux_rg[i, j] <- aux_rg[j, i] <- z[["rg"]]
    aux_se[i, j] <- aux_se[j, i] <- z[["se"]]
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(left = aux[[i]], right = aux[[j]], rg = z[["rg"]], se = z[["se"]])
  }
}

write.table(aux_rg, file.path(out_dir, "neurological_auxiliary_rg_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(aux_se, file.path(out_dir, "neurological_auxiliary_rg_se_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(do.call(rbind, pair_rows), file.path(out_dir, "neurological_auxiliary_rg_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gate_rows <- list()
target_rows <- list()
for (target in targets) {
  traits <- c(target, aux)
  mat <- diag(length(traits))
  dimnames(mat) <- list(traits, traits)
  mat[-1L, -1L] <- aux_rg
  for (j in seq_along(aux)) {
    z <- read_rg(target, aux[[j]])
    mat[1L, j + 1L] <- mat[j + 1L, 1L] <- z[["rg"]]
    target_rows[[length(target_rows) + 1L]] <- data.frame(target = target, auxiliary = aux[[j]], rg = z[["rg"]], se = z[["se"]])
  }
  aux_eig <- eigen(aux_rg, symmetric = TRUE, only.values = TRUE)$values
  full_eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  aux_cond <- max(aux_eig) / min(aux_eig)
  full_cond <- if (min(full_eig) > 0) max(full_eig) / min(full_eig) else Inf
  status <- if (min(aux_eig) <= 0 || aux_cond > 1e4 || min(full_eig) <= 0) "STOP" else if (aux_cond > 1e3 || full_cond > 1e3) "CONDITIONAL" else "PASS"
  gate_rows[[length(gate_rows) + 1L]] <- data.frame(
    target = target, n_aux = length(aux), aux_min_eigenvalue = min(aux_eig),
    aux_condition_number = aux_cond, full_min_eigenvalue = min(full_eig),
    full_condition_number = full_cond, matrix_gate = status
  )
  write.table(mat, file.path(out_dir, paste0(target, "_fixed_rg_matrix.tsv")), sep = "\t", quote = FALSE, col.names = NA)
}

write.table(do.call(rbind, target_rows), file.path(out_dir, "target_auxiliary_rg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, gate_rows), file.path(out_dir, "neurological_matrix_gate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("P5_NEURO_MATRIX_GATE_COMPLETE\n")

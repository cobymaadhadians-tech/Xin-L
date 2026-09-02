#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("target key required")
target_key <- args[[1L]]
targets <- c(
  "AD_GCST90704646_MAIN", "AD_GCST90704647_NOPROXY", "AD_GCST90704648_NOBIOBANK",
  "EPILEPSY_ILAE2023_EUR", "EPILEPSY_EHR_META_2OF3", "EPILEPSY_EHR_META_3OF3"
)
if (!target_key %in% targets) stop("unsupported target: ", target_key)

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})

aux <- c(
  "NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641",
  "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387"
)
traits <- c(target_key, aux)
paths <- c(
  setNames(file.path(root, "results/p5/ldsc_input", paste0(target_key, ".sumstats.gz")), target_key),
  setNames(file.path(root, "results/p5/neurological_auxiliary/ldsc_input", paste0(aux, ".sumstats.gz")), aux)
)
stopifnot(identical(names(paths), traits), all(file.exists(paths)))
munged <- lapply(paths, function(path) fread(path, showProgress = FALSE))
names(munged) <- traits

all_rg <- Cal_rg_h2g_alltraits(
  phenotype = traits,
  munged_sumstats = munged,
  ld_path = Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr")),
  wld_path = Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC")),
  sample_prev = NULL,
  population_prev = NULL
)
out_dir <- file.path(root, "results/p5/neurological_phbc/full_matrix", target_key)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(all_rg, file.path(out_dir, "ldsc_alltraits_observed.rds"))
write.table(all_rg$rg, file.path(out_dir, "ldsc_rg_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(all_rg$h2, file.path(out_dir, "ldsc_h2.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cat("P5_NEURO_FULL_MATRIX_COMPLETE\t", target_key, "\n", sep = "")

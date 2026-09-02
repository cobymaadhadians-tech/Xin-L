#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({library(data.table); library(pleioh2g)})
target <- "EPILEPSY_FINNGEN_R13"
aux <- c("NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641",
         "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387")
traits <- c(target, aux)
paths <- c(setNames(file.path(root, "results/p5/ldsc_input", paste0(target, ".sumstats.gz")), target),
           setNames(file.path(root, "results/p5/neurological_auxiliary/ldsc_input", paste0(aux, ".sumstats.gz")), aux))
if (!all(file.exists(paths))) stop("missing FinnGen or neurological auxiliary input")
munged <- lapply(paths, fread, showProgress = FALSE)
names(munged) <- traits
all_rg <- Cal_rg_h2g_alltraits(
  phenotype = traits, munged_sumstats = munged,
  ld_path = Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr")),
  wld_path = Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC")),
  sample_prev = NULL, population_prev = NULL
)
out <- file.path(root, "results/p6/epilepsy_finngen/full_matrix")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(all_rg, file.path(out, "ldsc_alltraits_observed.rds"))
write.table(all_rg$rg, file.path(out, "ldsc_rg_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(all_rg$h2, file.path(out, "ldsc_h2.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cat("P6_EPILEPSY_FINNGEN_FULL_MATRIX_COMPLETE\n")

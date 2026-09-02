#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

target <- "SCZ_PGC2022"
aux <- p25_panels[[target]]
traits <- c(target, aux)
paths <- p25_sumstats_paths(root)[traits]
paths[[target]] <- file.path(root, "results/p11_scz_neff_audit/input/SCZ_PGC2022_Nx2.sumstats.gz")
if (anyNA(paths) || !all(file.exists(paths))) stop("missing corrected matrix input")
munged <- lapply(paths, fread, showProgress = FALSE)
names(munged) <- traits

res <- Cal_rg_h2g_alltraits(
  phenotype = traits,
  munged_sumstats = munged,
  ld_path = Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr")),
  wld_path = Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC")),
  sample_prev = NULL,
  population_prev = NULL
)
out <- file.path(root, "results/p11_scz_neff_audit/full_matrix")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out, "ldsc_alltraits_observed_Nx2.rds"))
write.table(res$rg, file.path(out, "ldsc_rg_matrix_Nx2.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(res$gcov, file.path(out, "ldsc_gcov_matrix_Nx2.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(res$h2, file.path(out, "ldsc_h2_Nx2.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(res$h2Z, file.path(out, "ldsc_h2z_Nx2.tsv"), sep = "\t", quote = FALSE, col.names = NA)

old <- readRDS(file.path(root, "results/p25/full_matrix/SCZ_PGC2022/ldsc_alltraits_observed.rds"))
cmp <- data.table(
  metric = c("target_h2", "target_h2Z", "max_abs_rg_change", "target_gcov_scale", "max_aux_rg_change"),
  old = c(old$h2[1, target], old$h2Z[1, target], 0, 1, 0),
  corrected = c(
    res$h2[1, target], res$h2Z[1, target], max(abs(res$rg - old$rg)),
    median(res$gcov[target, ] / old$gcov[target, ]),
    max(abs(res$rg[aux, aux] - old$rg[aux, aux]))
  )
)
fwrite(cmp, file.path(out, "old_vs_Nx2_full_matrix.tsv"), sep = "\t")

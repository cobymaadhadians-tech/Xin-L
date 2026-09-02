#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
target <- args[[1]]
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))
allowed <- c("MDD_CLIN_PGC2025", "MDD_EHR_PGC2025", "MDD_QUEST_PGC2025", "MDD_FINNGEN_R13", "BD_CLIN_PGC4", "BD_FINNGEN_R13")
if (!target %in% allowed) stop("unsupported target")
aux <- p25_panels[[target]]
if (!"SCZ_PGC2022" %in% aux) stop("target panel does not contain PGC SCZ")
traits <- c(target, aux)
paths <- p25_sumstats_paths(root)[traits]
paths[["SCZ_PGC2022"]] <- file.path(root, "results/p11_scz_neff_audit/input/SCZ_PGC2022_Nx2.sumstats.gz")
munged <- lapply(paths, fread, showProgress = FALSE)
names(munged) <- traits
res <- Cal_rg_h2g_alltraits(
  phenotype = traits, munged_sumstats = munged,
  ld_path = Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr")),
  wld_path = Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC")),
  sample_prev = NULL, population_prev = NULL
)
out <- file.path(root, "results/p11_scz_neff_audit/aux_full_matrix", target)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out, "ldsc_alltraits_observed_Nx2.rds"))
old <- readRDS(file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds"))
cmp <- data.table(
  target = target,
  max_abs_rg_change = max(abs(res$rg - old$rg)),
  target_h2_old = old$h2[1, target], target_h2_corrected = res$h2[1, target],
  target_h2Z_old = old$h2Z[1, target], target_h2Z_corrected = res$h2Z[1, target],
  target_SCZ_rg_old = old$rg[target, "SCZ_PGC2022"],
  target_SCZ_rg_corrected = res$rg[target, "SCZ_PGC2022"]
)
fwrite(cmp, file.path(out, "old_vs_corrected.tsv"), sep = "\t")

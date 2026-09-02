#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out <- file.path(root, "results", "p11_genomicsem_scz_swap", "candidate_models")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

models <- list(
  broad_two_factor = paste(
    "F_INT =~ MDD + PTSD + ANX + OCD + AN",
    "F_EXT =~ SCZ + BD + ADHD + ASD + AUD + CUD",
    "F_INT ~~ F_EXT", sep = "\n"
  ),
  broad_two_factor_asd_int = paste(
    "F_INT =~ MDD + PTSD + ANX + OCD + AN + ASD",
    "F_EXT =~ SCZ + BD + ADHD + AUD + CUD",
    "F_INT ~~ F_EXT", sep = "\n"
  ),
  three_factor = paste(
    "F_INT =~ MDD + PTSD + ANX + OCD + AN",
    "F_PSY =~ SCZ + BD + ASD",
    "F_SUB =~ ADHD + AUD + CUD",
    "F_INT ~~ F_PSY + F_SUB",
    "F_PSY ~~ F_SUB", sep = "\n"
  )
)

checks <- list()
params <- list()
fits <- list()
for (target in c("SCZ_PGC2022", "SCZ_FINNGEN_R13")) {
  x <- readRDS(file.path(root, "results", "p11_genomicsem_scz_swap",
                         paste0(target, "_general_covstruc.rds")))
  for (id in names(models)) {
    cat("\nRUN", target, id, "\n")
    ans <- usermodel(
      covstruc = x$covstruc, estimation = "DWLS", model = models[[id]],
      CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE
    )
    ok <- !is.null(ans) && !is.null(ans$results) && nrow(as.data.frame(ans$results)) > 0L
    checks[[length(checks) + 1L]] <- data.frame(target_source = target, model_id = id, proper_solution = ok)
    if (ok) {
      p <- as.data.frame(ans$results)
      p$target_source <- target
      p$model_id <- id
      params[[length(params) + 1L]] <- p
      f <- as.data.frame(ans$modelfit)
      f$target_source <- target
      f$model_id <- id
      fits[[length(fits) + 1L]] <- f
      saveRDS(ans, file.path(out, paste0(target, "_", id, ".rds")))
    }
  }
}
fwrite(rbindlist(checks), file.path(out, "checks.tsv"), sep = "\t")
fwrite(rbindlist(params, fill = TRUE), file.path(out, "parameters.tsv"), sep = "\t")
fwrite(rbindlist(fits, fill = TRUE), file.path(out, "fit.tsv"), sep = "\t")
for (id in names(models)) writeLines(models[[id]], file.path(out, paste0(id, ".txt")))

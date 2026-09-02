#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[[1]])
if (!is.finite(task) || task < 0L || task > 19L) stop("task must be 0..19")

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out <- file.path(root, "results", "p11_genomicsem_scz_swap")
targets <- c("SCZ_PGC2022", "SCZ_FINNGEN_R13")
blocks <- (task * 10L + 1L):((task + 1L) * 10L)
model <- paste(
  "F_PSY =~ SCZ + BD + ASD",
  "F_SUB =~ ADHD + AUD + CUD",
  "F_PSY ~~ F_SUB",
  sep = "\n"
)

rows <- list()
fits <- list()
for (target in targets) {
  obj <- readRDS(file.path(out, paste0(target, "_primary_covstruc.rds")))
  for (b in blocks) {
    covstruc <- list(V = obj$covstruc$V, S = obj$blocks[, , b])
    fit <- tryCatch(
      usermodel(
        covstruc = covstruc, estimation = "DWLS", model = model,
        CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE
      ),
      error = function(e) e
    )
    if (inherits(fit, "error") || is.null(fit) || is.null(fit$results) ||
        nrow(as.data.frame(fit$results)) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        target_source = target, block = b, converged = FALSE,
        lhs = NA_character_, op = NA_character_, rhs = NA_character_,
        STD_All = NA_real_, error = if (inherits(fit, "error")) conditionMessage(fit) else "No interpretable parameter table"
      )
    } else {
      p <- as.data.frame(fit$results)
      p$target_source <- target
      p$block <- b
      p$converged <- TRUE
      p$error <- ""
      rows[[length(rows) + 1L]] <- p
      mf <- as.data.frame(fit$modelfit)
      mf$target_source <- target
      mf$block <- b
      fits[[length(fits) + 1L]] <- mf
    }
  }
}

fwrite(rbindlist(rows, fill = TRUE), file.path(out, sprintf("block_parameters_%02d.tsv", task)), sep = "\t")
if (length(fits)) {
  fwrite(rbindlist(fits, fill = TRUE), file.path(out, sprintf("block_fit_%02d.tsv", task)), sep = "\t")
} else {
  fwrite(data.table(), file.path(out, sprintf("block_fit_%02d.tsv", task)), sep = "\t")
}

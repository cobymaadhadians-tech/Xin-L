#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out <- file.path(root, "results", "p11_scz_neff_audit")
dir.create(file.path(out, "input"), recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  standardized = file.path(root, "results/p1/standardized_input/SCZ_PGC2022.standardized.tsv.gz"),
  ldsc = file.path(root, "results/p1/ldsc_input/SCZ_PGC2022.sumstats.gz")
)
outputs <- c(
  standardized = file.path(out, "input/SCZ_PGC2022_Nx2.standardized.tsv.gz"),
  ldsc = file.path(out, "input/SCZ_PGC2022_Nx2.sumstats.gz")
)

audit <- list()
for (kind in names(inputs)) {
  x <- fread(inputs[[kind]], showProgress = FALSE)
  if (!"N" %in% names(x) || any(!is.finite(x$N)) || any(x$N <= 0)) stop(kind, ": invalid source N")
  original <- x$N
  x[, N := 2 * N]
  if (!isTRUE(all.equal(x$N, 2 * original, tolerance = 0))) stop(kind, ": Nx2 validation failed")
  fwrite(x, outputs[[kind]], sep = "\t", compress = "gzip", quote = FALSE)
  audit[[kind]] <- data.table(
    input_kind = kind,
    rows = nrow(x),
    original_min = min(original), original_median = median(original), original_max = max(original),
    corrected_min = min(x$N), corrected_median = median(x$N), corrected_max = max(x$N),
    multiplier = 2,
    output = outputs[[kind]]
  )
}
fwrite(rbindlist(audit), file.path(out, "neff_input_audit.tsv"), sep = "\t")
writeLines(c(
  "source_field=PGC3_SCZ_NEFF",
  "source_header_description=effective sample size total",
  "source_median=58749.13",
  "analysis_convention=2_times_source_NEFF",
  "corrected_median=117498.26",
  "reason=PGC field follows half-balanced-effective-N convention for GenomicSEM/LDSC use"
), file.path(out, "neff_convention.txt"))

#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

n_block <- 200L
targets <- c("MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025")
out_dir <- file.path(root, "results/p9_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_target <- function(target) {
  traits <- c(target, p25_panels[[target]])
  full <- readRDS(file.path(root, "results/p25/full_matrix", target, "ldsc_alltraits_observed.rds"))
  files <- sort(list.files(
    file.path(root, "results/p25/phbc_cache", target),
    pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE
  ))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target)
  chunks <- lapply(files, readRDS)
  block_ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(block_ids) || !identical(sort(block_ids), seq_len(n_block))) {
    stop("incomplete or duplicated blocks for ", target)
  }
  rgarray <- array(
    NA_real_, c(length(traits), length(traits), n_block),
    dimnames = list(traits, traits, paste0("block_", seq_len(n_block)))
  )
  for (chunk in chunks) {
    if (!identical(chunk$traits, traits)) stop("trait-order mismatch for ", target)
    rgarray[, , chunk$block_ids] <- chunk$rgarray
  }
  list(target = target, traits = traits, full_rg = full$rg[traits, traits], rgarray = rgarray)
}

quadratic_form <- function(rg) {
  r <- rg[1L, -1L]
  aux <- rg[-1L, -1L, drop = FALSE]
  as.numeric(crossprod(r, solve(aux, r)))
}

inputs <- setNames(lapply(targets, load_target), targets)
records <- lapply(inputs, function(x) {
  list(
    point = quadratic_form(x$full_rg),
    blocks = vapply(seq_len(n_block), function(b) quadratic_form(x$rgarray[, , b]), numeric(1L))
  )
})

pairs <- data.frame(
  comparison = c("MDD_EHR_minus_QUEST", "MDD_FINNGEN_minus_QUEST"),
  left = c("MDD_EHR_PGC2025", "MDD_FINNGEN_R13"),
  right = c("MDD_QUEST_PGC2025", "MDD_QUEST_PGC2025"),
  stringsAsFactors = FALSE
)

rows <- lapply(seq_len(nrow(pairs)), function(i) {
  left <- records[[pairs$left[[i]]]]
  right <- records[[pairs$right[[i]]]]
  point <- left$point - right$point
  block_delta <- left$blocks - right$blocks
  pseudo <- n_block * point - (n_block - 1L) * block_delta
  se <- stats::sd(pseudo) / sqrt(n_block)
  z <- point / se
  data.frame(
    comparison = pairs$comparison[[i]],
    n_block = n_block,
    left_uncorrected_phbc = left$point,
    right_uncorrected_phbc = right$point,
    delta_fraction = point,
    delta_pp = 100 * point,
    paired_jackknife_se_fraction = se,
    paired_jackknife_se_pp = 100 * se,
    ci95_low_pp = 100 * (point - 1.96 * se),
    ci95_high_pp = 100 * (point + 1.96 * se),
    z = z,
    p_two_sided = 2 * stats::pnorm(-abs(z)),
    method = "correction-off 200-block paired delete-one-block jackknife",
    stringsAsFactors = FALSE
  )
})

write.table(
  do.call(rbind, rows),
  file.path(out_dir, "mdd_correction_off_paired_phbc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("P9_MDD_CORRECTION_OFF_COMPLETE\n")

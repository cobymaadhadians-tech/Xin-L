#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
base_dir <- file.path(root, "results/p8/p25_reestimated_weight_jackknife")
out_dir <- file.path(root, "results/p8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

targets <- c(
  "MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "MDD_QUEST_PGC2025",
  "SCZ_PGC2022", "SCZ_FINNGEN_R13"
)

read_target <- function(target) {
  files <- sort(list.files(
    file.path(base_dir, target), pattern = "^blocks_[0-9]{2}\\.tsv$", full.names = TRUE
  ))
  if (length(files) != 10L) stop("expected 10 chunks for ", target)
  out <- do.call(rbind, lapply(files, function(path) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  out <- out[order(out$block), ]
  if (!identical(out$block, seq_len(200L)) || anyDuplicated(out$block)) {
    stop("incomplete or duplicated blocks for ", target)
  }
  out
}

block_results <- do.call(rbind, lapply(targets, read_target))
write.table(
  block_results,
  file.path(out_dir, "p8_reestimated_weight_block_estimates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pairs <- data.frame(
  comparison = c(
    "MDD_EHR_minus_QUEST", "MDD_FINNGEN_minus_QUEST",
    "SCZ_FINNGEN_minus_PGC"
  ),
  left = c("MDD_EHR_PGC2025", "MDD_FINNGEN_R13", "SCZ_FINNGEN_R13"),
  right = c("MDD_QUEST_PGC2025", "MDD_QUEST_PGC2025", "SCZ_PGC2022"),
  stringsAsFactors = FALSE
)

rows <- vector("list", nrow(pairs))
records <- readRDS(file.path(root, "results/p25/final/p25_phbc_records.rds"))
for (i in seq_len(nrow(pairs))) {
  left <- block_results[block_results$target == pairs$left[[i]], ]
  right <- block_results[block_results$target == pairs$right[[i]], ]
  stopifnot(identical(left$block, right$block))

  delta_point <- as.numeric(records$full[[pairs$left[[i]]]]$phbc) -
    as.numeric(records$full[[pairs$right[[i]]]]$phbc)
  delta_jk <- left$reestimated_phbc - right$reestimated_phbc
  pseudo <- 200 * delta_point - 199 * delta_jk
  se <- stats::sd(pseudo) / sqrt(200)
  z <- delta_point / se
  fixed_delta_jk <- left$fixed_weight_phbc - right$fixed_weight_phbc
  fixed_pseudo <- 200 * delta_point - 199 * fixed_delta_jk
  fixed_se <- stats::sd(fixed_pseudo) / sqrt(200)

  rows[[i]] <- data.frame(
    comparison = pairs$comparison[[i]],
    n_block = 200L,
    delta_fraction = delta_point,
    delta_pp = 100 * delta_point,
    reestimated_weight_se_fraction = se,
    reestimated_weight_se_pp = 100 * se,
    ci95_low_pp = 100 * (delta_point - 1.96 * se),
    ci95_high_pp = 100 * (delta_point + 1.96 * se),
    z = z,
    p_two_sided = 2 * stats::pnorm(-abs(z)),
    fixed_weight_se_pp = 100 * fixed_se,
    se_change_pp = 100 * (se - fixed_se),
    se_ratio_reestimated_to_fixed = se / fixed_se,
    left_weight_min = min(left$reestimated_weight),
    left_weight_max = max(left$reestimated_weight),
    right_weight_min = min(right$reestimated_weight),
    right_weight_max = max(right$reestimated_weight),
    method = paste(
      "200-block paired delete-one-block jackknife; correction weight recalibrated",
      "at each leave-one-block-out point using the common 200-block covariance estimate"
    ),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, rows)
write.table(
  summary,
  file.path(out_dir, "p8_reestimated_weight_paired_delta_phbc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

weight_summary <- aggregate(
  reestimated_weight ~ target,
  data = block_results,
  FUN = function(v) c(min = min(v), median = stats::median(v), max = max(v), unique = length(unique(v)))
)
expanded <- data.frame(
  target = weight_summary$target,
  weight_min = weight_summary$reestimated_weight[, "min"],
  weight_median = weight_summary$reestimated_weight[, "median"],
  weight_max = weight_summary$reestimated_weight[, "max"],
  n_unique_weights = weight_summary$reestimated_weight[, "unique"],
  stringsAsFactors = FALSE
)
write.table(
  expanded,
  file.path(out_dir, "p8_reestimated_weight_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("P8_REESTIMATED_WEIGHT_COMBINE_COMPLETE\n")

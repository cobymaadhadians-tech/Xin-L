#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

n_block <- 200L
out_dir <- file.path(root, "results/p9_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_target <- function(target) {
  traits <- c(target, p25_panels[[target]])
  files <- sort(list.files(
    file.path(root, "results/p25/phbc_cache", target),
    pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE
  ))
  if (length(files) != 20L) stop("expected 20 cache chunks for ", target)
  chunks <- lapply(files, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_block))) {
    stop("incomplete or duplicated blocks for ", target)
  }
  rgarray <- array(
    NA_real_, c(length(traits), length(traits), n_block),
    dimnames = list(traits, traits, paste0("block_", seq_len(n_block)))
  )
  versions <- character(length(chunks))
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    if (!identical(chunk$traits, traits)) stop("trait-order mismatch for ", target)
    rgarray[, , chunk$block_ids] <- chunk$rgarray
    versions[[i]] <- chunk$pleioh2g_version
  }
  if (!all(is.finite(rgarray))) stop("non-finite rg cache for ", target)
  list(
    target = target,
    traits = traits,
    auxiliaries = p25_panels[[target]],
    rgarray = rgarray,
    versions = unique(versions),
    block_ids = ids
  )
}

lower_vector <- function(mat) mat[lower.tri(mat, diag = FALSE)]
inputs <- setNames(lapply(names(p25_panels), load_target), names(p25_panels))
hmp3 <- Sys.getenv("HM3_SNPLIST", unset = file.path(root, "reference/ldsc/w_hm3.snplist"))
hmp3_sha256 <- strsplit(system2("sha256sum", hmp3, stdout = TRUE), "[[:space:]]+")[[1L]][[1L]]
n_hmp3 <- length(readLines(hmp3, warn = FALSE)) - 1L

audit_rows <- lapply(inputs, function(x) {
  q <- length(x$traits) * (length(x$traits) - 1L) / 2L
  values <- vapply(seq_len(n_block), function(b) lower_vector(x$rgarray[, , b]), numeric(q))
  joint_covariance <- stats::cov(t(values)) * n_block
  eigenvalues <- eigen(joint_covariance, symmetric = TRUE, only.values = TRUE)$values
  offdiag <- joint_covariance[lower.tri(joint_covariance, diag = FALSE)]
  data.frame(
    target = x$target,
    n_traits = length(x$traits),
    n_aux = length(x$auxiliaries),
    n_blocks = n_block,
    block_ids_complete = identical(sort(x$block_ids), seq_len(n_block)),
    trait_order_consistent = TRUE,
    all_rg_finite = all(is.finite(x$rgarray)),
    all_rg_symmetric = all(vapply(seq_len(n_block), function(b) {
      isTRUE(all.equal(x$rgarray[, , b], t(x$rgarray[, , b]), tolerance = 1e-12))
    }, logical(1L))),
    all_rg_diagonal_one = all(vapply(seq_len(n_block), function(b) {
      max(abs(diag(x$rgarray[, , b]) - 1)) < 1e-12
    }, logical(1L))),
    joint_covariance_rows = nrow(joint_covariance),
    joint_covariance_cols = ncol(joint_covariance),
    expected_joint_dimension = q,
    joint_covariance_rank = qr(joint_covariance)$rank,
    joint_covariance_min_eigenvalue = min(eigenvalues),
    joint_covariance_max_eigenvalue = max(eigenvalues),
    offdiagonal_nonzero_fraction = mean(abs(offdiag) > 1e-15),
    max_abs_offdiagonal_covariance = max(abs(offdiag)),
    pleioh2g_version = paste(x$versions, collapse = ","),
    hmp3_snp_count = n_hmp3,
    hmp3_sha256 = hmp3_sha256,
    stringsAsFactors = FALSE
  )
})

families <- list(
  MDD = c("MDD_CLIN_PGC2025", "MDD_EHR_PGC2025", "MDD_QUEST_PGC2025", "MDD_FINNGEN_R13"),
  BD = c("BD_CLIN_PGC4", "BD_FINNGEN_R13"),
  SCZ = c("SCZ_PGC2022", "SCZ_FINNGEN_R13")
)
alignment_rows <- list()
for (family in names(families)) {
  combinations <- utils::combn(families[[family]], 2L, simplify = FALSE)
  for (pair in combinations) {
    left <- inputs[[pair[[1L]]]]
    right <- inputs[[pair[[2L]]]]
    if (!identical(left$auxiliaries, right$auxiliaries)) stop("family panel mismatch")
    aux <- left$auxiliaries
    difference <- left$rgarray[aux, aux, , drop = FALSE] - right$rgarray[aux, aux, , drop = FALSE]
    alignment_rows[[length(alignment_rows) + 1L]] <- data.frame(
      family = family,
      left_target = left$target,
      right_target = right$target,
      n_aux = length(aux),
      n_blocks = n_block,
      identical_block_ids = identical(sort(left$block_ids), sort(right$block_ids)),
      max_abs_auxiliary_rg_block_difference = max(abs(difference)),
      exactly_identical_auxiliary_rg_blocks = identical(
        left$rgarray[aux, aux, , drop = FALSE],
        right$rgarray[aux, aux, , drop = FALSE]
      ),
      stringsAsFactors = FALSE
    )
  }
}

write.table(
  do.call(rbind, audit_rows), file.path(out_dir, "joint_covariance_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, alignment_rows), file.path(out_dir, "auxiliary_block_alignment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("P9_JOINT_COVARIANCE_AUDIT_COMPLETE\n")

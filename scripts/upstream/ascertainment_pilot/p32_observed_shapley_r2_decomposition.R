#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 15)
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out_dir <- Sys.getenv("RESULTS_DIR", unset = file.path(root, "results", "p32_observed_shapley_r2"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
B <- 200L

models <- list(
  primary_six = c("BD_PGC2021", "MDD_CLIN_PGC2025", "ADHD", "OCD_2025", "PTSD", "AUD"),
  sensitivity_all_ten = c(
    "MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD",
    "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"
  )
)

source_specs <- list(
  PGC = list(
    target = "SCZ_PGC2022",
    full = file.path(root, "results", "p11_scz_neff_audit", "full_matrix",
                      "ldsc_alltraits_observed_Nx2.rds"),
    chunks = file.path(root, "results", "p11_scz_neff_audit", "phbc_cache")
  ),
  narrow = list(
    target = "SCZ_FINNGEN_R13",
    full = file.path(root, "results", "p25", "full_matrix", "SCZ_FINNGEN_R13",
                     "ldsc_alltraits_observed.rds"),
    chunks = file.path(root, "results", "p25", "phbc_cache", "SCZ_FINNGEN_R13")
  ),
  broad = list(
    target = "SCZ_FINNGEN_R13_EXMORE",
    full = file.path(root, "results", "p26_broad_scz", "full_matrix",
                     "SCZ_FINNGEN_R13_EXMORE", "ldsc_alltraits_observed.rds"),
    chunks = file.path(root, "results", "p26_broad_scz", "phbc_cache",
                       "SCZ_FINNGEN_R13_EXMORE")
  )
)

bits <- function(mask, p) {
  if (mask == 0L) return(integer())
  which(as.logical(intToBits(mask)[seq_len(p)]))
}

read_source <- function(spec, source_name) {
  if (!file.exists(spec$full)) stop("Missing full matrix: ", spec$full)
  paths <- sort(list.files(spec$chunks, "^chunk_[0-9]{3}\\.rds$", full.names = TRUE))
  if (length(paths) != 20L) stop(source_name, ": expected 20 chunks, found ", length(paths))
  full <- readRDS(spec$full)
  R <- full$rg
  if (is.null(R) || !identical(rownames(R), colnames(R))) {
    stop(source_name, ": invalid full rg matrix")
  }
  chunks <- lapply(paths, readRDS)
  ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (!identical(sort(as.integer(ids)), seq_len(B))) {
    stop(source_name, ": block coverage is incomplete or duplicated")
  }
  n <- nrow(R)
  blocks <- array(NA_real_, c(n, n, B),
                  dimnames = list(rownames(R), colnames(R), paste0("block_", seq_len(B))))
  for (ch in chunks) {
    if (is.null(ch$rgarray) || !identical(dim(ch$rgarray)[1:2], c(n, n))) {
      stop(source_name, ": invalid block rgarray")
    }
    blocks[, , ch$block_ids] <- ch$rgarray
  }
  if (any(!is.finite(R)) || any(!is.finite(blocks))) {
    stop(source_name, ": non-finite rg values")
  }
  list(target = spec$target, full = R, blocks = blocks)
}

subset_r2 <- function(R, target, predictors) {
  p <- length(predictors)
  values <- numeric(2^p)
  for (mask in 1:(2^p - 1L)) {
    members <- predictors[bits(mask, p)]
    idx <- c(target, members)
    Rt <- R[idx, idx, drop = FALSE]
    Rx <- Rt[-1L, -1L, drop = FALSE]
    rxy <- Rt[-1L, 1L, drop = FALSE]
    if (min(eigen(Rx, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop("Non-positive-definite predictor matrix at mask ", mask)
    }
    values[mask + 1L] <- as.numeric(crossprod(rxy, solve(Rx, rxy)))
  }
  if (any(!is.finite(values)) || any(values < -1e-8) || any(values > 1 + 1e-8)) {
    stop("Invalid subset R2")
  }
  pmax(0, pmin(1, values))
}

shapley <- function(values, predictors) {
  p <- length(predictors)
  phi <- numeric(p)
  names(phi) <- predictors
  for (j in seq_len(p)) {
    others <- setdiff(seq_len(p), j)
    for (submask in 0:(2^(p - 1L) - 1L)) {
      members <- others[bits(submask, p - 1L)]
      without <- if (length(members)) sum(2^(members - 1L)) else 0L
      with <- without + 2^(j - 1L)
      k <- length(members)
      weight <- factorial(k) * factorial(p - k - 1L) / factorial(p)
      phi[j] <- phi[j] + weight * (values[with + 1L] - values[without + 1L])
    }
  }
  phi
}

compute_model <- function(source, source_name, model_name, predictors) {
  p <- length(predictors)
  full_values <- subset_r2(source$full, source$target, predictors)
  full_phi <- shapley(full_values, predictors)
  if (abs(sum(full_phi) - full_values[2^p]) > 1e-8) {
    stop(source_name, "/", model_name, ": full Shapley sum check failed")
  }
  block_phi <- matrix(NA_real_, B, p, dimnames = list(seq_len(B), predictors))
  block_r2 <- numeric(B)
  for (b in seq_len(B)) {
    vals <- subset_r2(source$blocks[, , b], source$target, predictors)
    ph <- shapley(vals, predictors)
    if (abs(sum(ph) - vals[2^p]) > 1e-8) {
      stop(source_name, "/", model_name, ": block Shapley sum check failed at block ", b)
    }
    block_phi[b, ] <- ph
    block_r2[b] <- vals[2^p]
    if (b %% 25L == 0L) cat(source_name, model_name, "block", b, "of", B, "\n")
  }
  list(full_values = full_values, full_r2 = full_values[2^p],
       full_phi = full_phi, block_r2 = block_r2, block_phi = block_phi)
}

jk <- function(full, blocks) {
  full <- unname(full)
  blocks <- unname(blocks)
  center <- mean(blocks)
  se <- sqrt((B - 1) / B * sum((blocks - center)^2))
  z <- full / se
  c(estimate = full, se = se, CI95_low = full - qnorm(.975) * se,
    CI95_high = full + qnorm(.975) * se, z = z,
    P = 2 * pnorm(-abs(z)), min_delete_one = min(blocks),
    max_delete_one = max(blocks),
    direction_reversals_200 = sum(sign(blocks) != sign(full)))
}

profile_test <- function(full_delta, block_delta) {
  p <- length(full_delta)
  keep <- seq_len(p - 1L)
  X <- sweep(block_delta[, keep, drop = FALSE], 2L,
             colMeans(block_delta[, keep, drop = FALSE]), "-")
  V <- (B - 1) / B * crossprod(X)
  stat <- as.numeric(crossprod(full_delta[keep], solve(V, full_delta[keep])))
  c(statistic = stat, df = p - 1L,
    P = pchisq(stat, df = p - 1L, lower.tail = FALSE),
    covariance_rank = qr(V)$rank)
}

profile_test_full <- function(full_delta, block_delta) {
  p <- length(full_delta)
  X <- sweep(block_delta, 2L, colMeans(block_delta), "-")
  V <- (B - 1) / B * crossprod(X)
  rk <- qr(V)$rank
  cond <- if (rk == p) kappa(V, exact = TRUE) else NA_real_
  if (rk < p) {
    return(c(statistic = NA_real_, df = NA_real_, P = NA_real_,
             covariance_rank = rk, condition_number = cond))
  }
  stat <- as.numeric(crossprod(full_delta, solve(V, full_delta)))
  c(statistic = stat, df = p,
    P = pchisq(stat, df = p, lower.tail = FALSE),
    covariance_rank = rk, condition_number = cond)
}

cat("Loading PGC, narrow FinnGen and broad FinnGen matrices...\n")
sources <- lapply(names(source_specs), function(nm) read_source(source_specs[[nm]], nm))
names(sources) <- names(source_specs)
all_aux <- unique(unlist(models, use.names = FALSE))
aux_validation <- do.call(rbind, lapply(names(sources)[-1L], function(nm) {
  data.frame(
    reference = "PGC", source = nm,
    full_max_abs_auxiliary_difference = max(abs(
      sources[[nm]]$full[all_aux, all_aux] - sources$PGC$full[all_aux, all_aux])),
    block_max_abs_auxiliary_difference = max(vapply(seq_len(B), function(b) max(abs(
      sources[[nm]]$blocks[all_aux, all_aux, b] -
        sources$PGC$blocks[all_aux, all_aux, b])), numeric(1)))
  )
}))
if (any(as.matrix(aux_validation[, -c(1L, 2L)]) > 1e-8)) {
  stop("Auxiliary-only matrices are not matched")
}

results <- list()
full_rows <- list()
subset_rows <- list()
block_rows <- list()
for (model_name in names(models)) {
  predictors <- models[[model_name]]
  p <- length(predictors)
  for (source_name in names(sources)) {
    cat("Computing", source_name, model_name, "\n")
    ans <- compute_model(sources[[source_name]], source_name, model_name, predictors)
    results[[paste(source_name, model_name, sep = "__")]] <- ans
    full_rows[[length(full_rows) + 1L]] <- data.frame(
      model = model_name, source = source_name, predictor = predictors,
      full_R2 = ans$full_r2, shapley_phi = ans$full_phi,
      shapley_share_of_R2 = ans$full_phi / ans$full_r2,
      sum_phi_minus_R2 = sum(ans$full_phi) - ans$full_r2
    )
    masks <- 0:(2^p - 1L)
    subset_rows[[length(subset_rows) + 1L]] <- data.frame(
      model = model_name, source = source_name, mask = masks,
      subset_size = vapply(masks, function(m) length(bits(m, p)), integer(1)),
      subset_predictors = vapply(masks, function(m) {
        z <- predictors[bits(m, p)]
        if (length(z)) paste(z, collapse = "+") else "(intercept only)"
      }, character(1)),
      R2 = unname(ans$full_values)
    )
    block_rows[[length(block_rows) + 1L]] <- do.call(rbind, lapply(seq_len(B), function(b) {
      data.frame(model = model_name, source = source_name, block = b,
                 predictor = predictors, R2 = ans$block_r2[b],
                 shapley_phi = ans$block_phi[b, ],
                 shapley_share = ans$block_phi[b, ] / ans$block_r2[b])
    }))
  }
}

full_table <- do.call(rbind, full_rows)
subset_table <- do.call(rbind, subset_rows)
block_table <- do.call(rbind, block_rows)
comparisons <- list(
  narrow_minus_PGC = c("narrow", "PGC"),
  broad_minus_PGC = c("broad", "PGC"),
  broad_minus_narrow = c("broad", "narrow")
)
profile_rows <- list()
omnibus_rows <- list()
alternative_omnibus_rows <- list()
total_rows <- list()
for (model_name in names(models)) {
  predictors <- models[[model_name]]
  for (comparison in names(comparisons)) {
    left <- comparisons[[comparison]][1L]
    right <- comparisons[[comparison]][2L]
    lk <- paste(left, model_name, sep = "__")
    rk <- paste(right, model_name, sep = "__")
    delta <- results[[lk]]$full_phi - results[[rk]]$full_phi
    delta_blocks <- results[[lk]]$block_phi - results[[rk]]$block_phi
    profile_rows[[length(profile_rows) + 1L]] <- do.call(rbind, lapply(seq_along(predictors), function(j) {
      inf <- jk(delta[j], delta_blocks[, j])
      data.frame(model = model_name, comparison = comparison,
                 predictor = predictors[j], phi_left = results[[lk]]$full_phi[j],
                 phi_right = results[[rk]]$full_phi[j],
                 delta_phi_left_minus_right = unname(inf["estimate"]),
                 paired_jackknife_SE = unname(inf["se"]),
                 CI95_low = unname(inf["CI95_low"]),
                 CI95_high = unname(inf["CI95_high"]),
                 z = unname(inf["z"]), P = unname(inf["P"]),
                 min_delete_one_delta = unname(inf["min_delete_one"]),
                 max_delete_one_delta = unname(inf["max_delete_one"]),
                 direction_reversals_200 = as.integer(inf["direction_reversals_200"]))
    }))
    om <- profile_test(delta, delta_blocks)
    omnibus_rows[[length(omnibus_rows) + 1L]] <- data.frame(
      model = model_name, comparison = comparison,
      statistic = unname(om["statistic"]), df = unname(om["df"]),
      P = unname(om["P"]), covariance_rank = unname(om["covariance_rank"]),
      omitted_component_for_rank_constraint = tail(predictors, 1L)
    )

    om_full <- profile_test_full(delta, delta_blocks)
    alternative_omnibus_rows[[length(alternative_omnibus_rows) + 1L]] <- data.frame(
      model = model_name, comparison = comparison,
      test = "raw_phi_full_rank",
      statistic = unname(om_full["statistic"]), df = unname(om_full["df"]),
      P = unname(om_full["P"]), covariance_rank = unname(om_full["covariance_rank"]),
      condition_number = unname(om_full["condition_number"]),
      omitted_component = NA_character_
    )

    left_phi_blocks <- results[[lk]]$block_phi
    right_phi_blocks <- results[[rk]]$block_phi
    left_share_blocks <- left_phi_blocks / results[[lk]]$block_r2
    right_share_blocks <- right_phi_blocks / results[[rk]]$block_r2
    delta_share <- results[[lk]]$full_phi / results[[lk]]$full_r2 -
      results[[rk]]$full_phi / results[[rk]]$full_r2
    delta_share_blocks <- left_share_blocks - right_share_blocks
    om_share <- profile_test(delta_share, delta_share_blocks)
    alternative_omnibus_rows[[length(alternative_omnibus_rows) + 1L]] <- data.frame(
      model = model_name, comparison = comparison,
      test = "normalized_share_rank_constrained",
      statistic = unname(om_share["statistic"]), df = unname(om_share["df"]),
      P = unname(om_share["P"]), covariance_rank = unname(om_share["covariance_rank"]),
      condition_number = NA_real_,
      omitted_component = tail(predictors, 1L)
    )

    total_inf <- jk(results[[lk]]$full_r2 - results[[rk]]$full_r2,
                    results[[lk]]$block_r2 - results[[rk]]$block_r2)
    total_rows[[length(total_rows) + 1L]] <- data.frame(
      model = model_name, comparison = comparison,
      R2_left = results[[lk]]$full_r2, R2_right = results[[rk]]$full_r2,
      delta_R2_left_minus_right = unname(total_inf["estimate"]),
      paired_jackknife_SE = unname(total_inf["se"]),
      CI95_low = unname(total_inf["CI95_low"]),
      CI95_high = unname(total_inf["CI95_high"]),
      z = unname(total_inf["z"]), P = unname(total_inf["P"]),
      min_delete_one_delta = unname(total_inf["min_delete_one"]),
      max_delete_one_delta = unname(total_inf["max_delete_one"]),
      direction_reversals_200 = as.integer(total_inf["direction_reversals_200"])
    )
  }
}
profile_table <- do.call(rbind, profile_rows)
profile_table$BH_q_within_comparison_model <- ave(
  profile_table$P, profile_table$comparison, profile_table$model,
  FUN = function(x) p.adjust(x, "BH")
)
omnibus_table <- do.call(rbind, omnibus_rows)
omnibus_table$BH_q_within_model <- ave(
  omnibus_table$P, omnibus_table$model, FUN = function(x) p.adjust(x, "BH")
)
alternative_omnibus_table <- do.call(rbind, alternative_omnibus_rows)
alternative_omnibus_table$BH_q_within_model_test <- ave(
  alternative_omnibus_table$P, alternative_omnibus_table$test,
  alternative_omnibus_table$model, FUN = function(x) p.adjust(x, "BH")
)
total_table <- do.call(rbind, total_rows)
total_table$BH_q_within_model <- ave(
  total_table$P, total_table$model, FUN = function(x) p.adjust(x, "BH")
)

write.table(full_table, file.path(out_dir, "shapley_full.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(subset_table, file.path(out_dir, "subset_R2_full.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(block_table, file.path(out_dir, "shapley_delete_one_blocks.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(profile_table, file.path(out_dir, "shapley_profile_contrasts.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(omnibus_table, file.path(out_dir, "shapley_profile_omnibus.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(alternative_omnibus_table,
            file.path(out_dir, "shapley_profile_omnibus_alternatives.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(total_table, file.path(out_dir, "shapley_total_R2_contrasts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aux_validation, file.path(out_dir, "matched_auxiliary_validation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(list(full = full_table, subsets = subset_table, blocks = block_table,
             profile = profile_table, omnibus = omnibus_table,
             alternative_omnibus = alternative_omnibus_table,
             total_R2 = total_table,
             auxiliary_validation = aux_validation),
        file.path(out_dir, "shapley_r2_complete.rds"))
writeLines(c(
  "analysis = Shapley/general-dominance decomposition of multivariable genetic R2",
  "sources = PGC, narrow FinnGen, broad FinnGen",
  "models = primary six predictors and sensitivity ten predictors",
  "uncertainty = 200 matched delete-one-block paired jackknife",
  "profile_omnibus = raw-phi rank-constrained, raw-phi full-rank and normalized-share tests",
  "total_R2_contrasts = paired delete-one-block jackknife",
  paste0("max_auxiliary_matrix_difference = ", format(max(as.matrix(aux_validation[, -c(1L, 2L)])), scientific = TRUE)),
  paste0("max_full_sum_validation_error = ", format(max(abs(full_table$sum_phi_minus_R2)), scientific = TRUE)),
  paste0("completed_utc = ", format(Sys.time(), tz = "UTC")),
  "status = COMPLETED"
), file.path(out_dir, "STATUS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("COMPLETED\n")
print(full_table[, c("model", "source", "predictor", "full_R2", "shapley_phi")])
print(profile_table[, c("model", "comparison", "predictor",
                        "delta_phi_left_minus_right", "paired_jackknife_SE", "P",
                        "BH_q_within_comparison_model")])
print(omnibus_table)
print(alternative_omnibus_table)
print(total_table)

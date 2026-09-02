#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out <- file.path(root, "results", "p11_genomicsem_scz_swap")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

targets <- c("SCZ_PGC2022", "SCZ_FINNGEN_R13")
name_map <- c(
  SCZ_PGC2022 = "SCZ", SCZ_FINNGEN_R13 = "SCZ",
  MDD_CLIN_PGC2025 = "MDD", BD_PGC2021 = "BD", PTSD = "PTSD",
  ANX = "ANX", ADHD = "ADHD", ASD = "ASD", OCD_2025 = "OCD",
  AN = "AN", AUD = "AUD", CUD_2023_EUR = "CUD"
)

primary_traits <- c("BD_PGC2021", "ASD", "ADHD", "AUD", "CUD_2023_EUR")
primary_model <- paste(
  "F_PSY =~ SCZ + BD + ASD",
  "F_SUB =~ ADHD + AUD + CUD",
  "F_PSY ~~ F_SUB",
  sep = "\n"
)

general_aux <- c(
  "MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD", "ASD",
  "OCD_2025", "AN", "AUD", "CUD_2023_EUR"
)
general_model <- "F_P =~ SCZ + MDD + BD + PTSD + ANX + ADHD + ASD + OCD + AN + AUD + CUD"

vech <- function(m) m[lower.tri(m, diag = TRUE)]

read_target <- function(target) {
  corrected_pgc <- identical(target, "SCZ_PGC2022")
  full_path <- if (corrected_pgc) {
    file.path(root, "results", "p11_scz_neff_audit", "full_matrix", "ldsc_alltraits_observed_Nx2.rds")
  } else {
    file.path(root, "results", "p25", "full_matrix", target, "ldsc_alltraits_observed.rds")
  }
  full <- readRDS(full_path)
  chunk_dir <- if (corrected_pgc) {
    file.path(root, "results", "p11_scz_neff_audit", "phbc_cache")
  } else {
    file.path(root, "results", "p25", "phbc_cache", target)
  }
  chunk_paths <- sort(list.files(
    chunk_dir,
    pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE
  ))
  if (length(chunk_paths) != 20L) stop(target, ": expected 20 cache chunks, found ", length(chunk_paths))
  chunks <- lapply(chunk_paths, readRDS)
  block_ids <- unlist(lapply(chunks, `[[`, "block_ids"), use.names = FALSE)
  if (!identical(sort(block_ids), 1:200)) stop(target, ": cache blocks are incomplete or duplicated")
  arrays <- lapply(chunks, `[[`, "gcovarray")
  all_names <- rownames(full$gcov)
  arr <- array(NA_real_, dim = c(length(all_names), length(all_names), 200L),
               dimnames = list(all_names, all_names, paste0("block_", 1:200)))
  for (i in seq_along(chunks)) arr[, , chunks[[i]]$block_ids] <- arrays[[i]]
  list(full = full$gcov, h2 = full$h2, h2Z = full$h2Z, blocks = arr)
}

make_covstruc <- function(dat, target, aux) {
  original <- c(target, aux)
  S <- dat$full[original, original, drop = FALSE]
  blocks <- dat$blocks[original, original, , drop = FALSE]
  short <- unname(name_map[original])
  if (anyNA(short)) stop("Missing short name for: ", paste(original[is.na(short)], collapse = ", "))
  dimnames(S) <- list(short, short)
  dimnames(blocks)[1:2] <- list(short, short)
  theta <- t(vapply(seq_len(dim(blocks)[3]), function(b) vech(blocks[, , b]), numeric(length(short) * (length(short) + 1) / 2)))
  centered <- sweep(theta, 2, colMeans(theta), "-")
  B <- nrow(theta)
  V <- (B - 1) / B * crossprod(centered)
  diag_positions <- cumsum(c(1L, rev(seq_len(length(short) - 1L)) + 1L))
  diag_se <- sqrt(diag(V))[diag_positions]
  observed_se <- unname(dat$h2[1, original] / dat$h2Z[1, original])
  validation <- data.frame(
    trait = short,
    jackknife_h2_se = diag_se,
    reported_h2_se = observed_se,
    ratio = diag_se / observed_se
  )
  # GenomicSEM::usermodel() reads covstruc positionally: [[1]] is V and
  # [[2]] is S. Preserve the native ldsc() object order.
  list(covstruc = list(V = V, S = S), blocks = blocks, validation = validation,
       V_min_eigenvalue = min(eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values),
       V_condition_number = kappa(V))
}

read_metric <- function(x, key) {
  y <- as.data.frame(x)
  if (key %in% names(y)) return(suppressWarnings(as.numeric(y[[key]][1])))
  NA_real_
}

fit_one <- function(covstruc, model, target, model_id) {
  fit <- usermodel(
    covstruc = covstruc,
    estimation = "DWLS",
    model = model,
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = TRUE,
    fix_resid = TRUE
  )
  if (is.null(fit) || is.null(fit$results) || nrow(as.data.frame(fit$results)) == 0L) {
    stop("GenomicSEM returned no interpretable parameter table")
  }
  params <- as.data.frame(fit$results)
  params$target_source <- target
  params$model_id <- model_id
  mf <- as.data.frame(fit$modelfit)
  mf$target_source <- target
  mf$model_id <- model_id
  list(fit = fit, params = params, fit_metrics = mf)
}

all_params <- list()
all_fit <- list()
checks <- list()

for (target in targets) {
  dat <- read_target(target)
  primary <- make_covstruc(dat, target, primary_traits)
  general <- make_covstruc(dat, target, general_aux)

  saveRDS(primary, file.path(out, paste0(target, "_primary_covstruc.rds")))
  saveRDS(general, file.path(out, paste0(target, "_general_covstruc.rds")))
  fwrite(primary$validation, file.path(out, paste0(target, "_primary_V_validation.tsv")), sep = "\t")
  fwrite(general$validation, file.path(out, paste0(target, "_general_V_validation.tsv")), sep = "\t")

  for (spec in list(
    list(id = "primary_two_factor", obj = primary, model = primary_model),
    list(id = "general_one_factor", obj = general, model = general_model)
  )) {
    ans <- tryCatch(
      fit_one(spec$obj$covstruc, spec$model, target, spec$id),
      error = function(e) e
    )
    if (inherits(ans, "error")) {
      checks[[length(checks) + 1L]] <- data.frame(
        target_source = target, model_id = spec$id, converged = FALSE,
        error = conditionMessage(ans), V_min_eigenvalue = spec$obj$V_min_eigenvalue,
        V_condition_number = spec$obj$V_condition_number
      )
    } else {
      saveRDS(ans$fit, file.path(out, paste0(target, "_", spec$id, "_fit.rds")))
      all_params[[length(all_params) + 1L]] <- ans$params
      all_fit[[length(all_fit) + 1L]] <- ans$fit_metrics
      checks[[length(checks) + 1L]] <- data.frame(
        target_source = target, model_id = spec$id, converged = TRUE,
        error = "", V_min_eigenvalue = spec$obj$V_min_eigenvalue,
        V_condition_number = spec$obj$V_condition_number
      )
    }
  }
}

fwrite(rbindlist(all_params, fill = TRUE), file.path(out, "full_model_parameters.tsv"), sep = "\t")
fwrite(rbindlist(all_fit, fill = TRUE), file.path(out, "full_model_fit.tsv"), sep = "\t")
fwrite(rbindlist(checks, fill = TRUE), file.path(out, "full_model_checks.tsv"), sep = "\t")
writeLines(primary_model, file.path(out, "primary_two_factor_model.txt"))
writeLines(general_model, file.path(out, "general_one_factor_model.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

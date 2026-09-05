suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("ANALYSIS_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), mustWork = FALSE)
out_dir <- file.path(root, "results/d6b_direct_contrast_20260905")
capture_path <- Sys.getenv("D6B_CAPTURE_RDS", file.path(out_dir, "D6B_direct_textinject_capture.rds"))
sem_path <- file.path(root, "results/d6b_published_pi_braincell_stratified_gsem_20260905/D6B_S_LDSCOutput.rds")
enrichment_path <- file.path(root, "results/d6b_published_pi_braincell_stratified_gsem_20260905/D6B_enrichment_long.tsv")
stopifnot(file.exists(capture_path), file.exists(sem_path), file.exists(enrichment_path))

cap <- readRDS(capture_path)
s <- readRDS(sem_path)
en <- fread(enrichment_path)

parameter_map <- c(
  "F_COMP~~SCZ_PGC", "F_INT~~SCZ_PGC", "F_SUD~~SCZ_PGC",
  "F_COMP~~SCZ_FG", "F_INT~~SCZ_FG", "F_SUD~~SCZ_FG"
)
custom <- unique(en[annotation_class == "custom",
  .(Annotation, annotation_key, annotation_basename, annotation_class, published_label)])
qc <- en[, .(
  Annotation_model_warning = any(Annotation_model_warning %in% TRUE | Warning %in% c(TRUE, "TRUE", "1", "1.0")),
  Annotation_model_error = any(Annotation_model_error %in% TRUE | Error %in% c(TRUE, "TRUE", "1", "1.0"))
), by = Annotation]
custom <- merge(custom, qc, by = "Annotation", all.x = TRUE, sort = FALSE)
setorder(custom, Annotation)
stopifnot(nrow(custom) == 29L)

match_capture <- function(annotation) {
  ii <- which(vapply(cap, function(x) {
    identical(as.character(names(s$S)[x$n]), annotation)
  }, logical(1)))
  if (length(ii)) cap[[ii[1L]]] else NULL
}

make_record <- function(annotation, i, cc, raw, covariance, status, z_smooth_pgc, z_smooth_finngen) {
  domains <- c("COMP", "INT", "SUD")
  delta <- raw[4:6] - raw[1:3]
  contrast_matrix <- cbind(-diag(3), diag(3))
  delta_covariance <- contrast_matrix %*% covariance %*% t(contrast_matrix)
  variance <- diag(delta_covariance)
  se <- ifelse(is.finite(variance) & variance >= 0, sqrt(variance), NA_real_)
  z <- ifelse(is.finite(delta) & is.finite(se) & se > 0, delta / se, NA_real_)
  p <- ifelse(is.finite(z), 2 * pnorm(-abs(z)), NA_real_)
  finite <- is.finite(delta) & is.finite(se) & is.finite(z) & is.finite(p)
  list(
    domain = domains[i],
    raw_covariance_PGC = raw[i],
    raw_covariance_FinnGen = raw[i + 3L],
    delta_FinnGen_minus_PGC = delta[i],
    SE = se[i],
    z = z[i],
    P_two_sided = p[i],
    contrast_status = status,
    finite_contrast = finite[i],
    Z_smooth_PGC = z_smooth_pgc,
    Z_smooth_FinnGen = z_smooth_finngen,
    max_abs_Z_smooth = if (all(is.finite(c(z_smooth_pgc, z_smooth_finngen)))) {
      max(abs(c(z_smooth_pgc, z_smooth_finngen)))
    } else NA_real_,
    covariance_condition_number = if (all(is.finite(delta_covariance))) {
      ev <- eigen(delta_covariance, symmetric = TRUE, only.values = TRUE)$values
      if (all(ev > 0)) max(ev) / min(ev) else NA_real_
    } else NA_real_,
    delta_covariance = delta_covariance
  )
}

raw_rows <- vector("list", nrow(custom) * 3L)
k <- 0L

for (j in seq_len(nrow(custom))) {
  annotation <- custom$Annotation[j]
  cc <- match_capture(annotation)
  status <- if (is.null(cc) || is.null(cc$covariance)) "no sandwich capture" else "sandwich captured"
  raw <- rep(NA_real_, 6L)
  covariance <- matrix(NA_real_, 6L, 6L)
  if (!is.null(cc) && !is.null(cc$covariance) && length(cc$estimate) == 6L) {
    idx <- match(parameter_map, cc$map)
    if (all(!is.na(idx))) {
      raw <- as.numeric(cc$estimate)
      covariance <- as.matrix(cc$covariance)[idx, idx, drop = FALSE]
    }
  }
  z_smooth <- vapply(parameter_map, function(parameter) {
    value <- en[Annotation == annotation & target_parameter == parameter, Z_smooth]
    if (length(value)) as.numeric(value[1L]) else NA_real_
  }, numeric(1))
  records <- lapply(seq_len(3L), function(i) {
    make_record(annotation, i, cc, raw, covariance, status,
                z_smooth_pgc = z_smooth[i],
                z_smooth_finngen = z_smooth[i + 3L])
  })
  for (i in seq_len(3L)) {
    k <- k + 1L
    r <- records[[i]]
    raw_rows[[k]] <- data.table(
      Annotation = annotation,
      annotation_key = custom$annotation_key[j],
      annotation_basename = custom$annotation_basename[j],
      annotation_class = custom$annotation_class[j],
      published_label = custom$published_label[j],
      domain = r$domain,
      raw_covariance_PGC = r$raw_covariance_PGC,
      raw_covariance_FinnGen = r$raw_covariance_FinnGen,
      delta_FinnGen_minus_PGC = r$delta_FinnGen_minus_PGC,
      SE = r$SE,
      z = r$z,
      P_two_sided = r$P_two_sided,
      alpha_domain_29 = 0.05 / nrow(custom),
      Bonferroni_significant_domain_29 = isTRUE(is.finite(r$P_two_sided) && r$P_two_sided <= 0.05 / nrow(custom)),
      QC_model_warning_free = isTRUE(r$finite_contrast) &&
        !isTRUE(custom$Annotation_model_warning[j]) &&
        !isTRUE(custom$Annotation_model_error[j]),
      QC_zsmooth_ok = isTRUE(r$finite_contrast) && is.finite(r$max_abs_Z_smooth) && r$max_abs_Z_smooth <= 1.96,
      QC_primary_eligible = isTRUE(r$finite_contrast) &&
        !isTRUE(custom$Annotation_model_warning[j]) &&
        !isTRUE(custom$Annotation_model_error[j]) &&
        is.finite(r$max_abs_Z_smooth) && r$max_abs_Z_smooth <= 1.96,
      Bonferroni_significant_domain_29_primary = isTRUE(
        is.finite(r$P_two_sided) && r$P_two_sided <= 0.05 / nrow(custom) &&
          r$finite_contrast && !isTRUE(custom$Annotation_model_warning[j]) &&
          !isTRUE(custom$Annotation_model_error[j]) &&
          is.finite(r$max_abs_Z_smooth) && r$max_abs_Z_smooth <= 1.96
      ),
      Annotation_model_warning = custom$Annotation_model_warning[j],
      Annotation_model_error = custom$Annotation_model_error[j],
      Z_smooth_PGC = r$Z_smooth_PGC,
      Z_smooth_FinnGen = r$Z_smooth_FinnGen,
      max_abs_Z_smooth = r$max_abs_Z_smooth,
      contrast_status = r$contrast_status,
      finite_contrast = r$finite_contrast,
      covariance_condition_number = r$covariance_condition_number
    )
  }

}

raw_table <- rbindlist(raw_rows, fill = TRUE)
fwrite(raw_table, file.path(out_dir, "D6B_raw_annotation_covariance_contrast.tsv"), sep = "\t", na = "NA", quote = FALSE)
domain_summary <- raw_table[, .(
  n_rows = .N,
  n_finite = sum(finite_contrast %in% TRUE),
  n_significant_domain_29 = sum(Bonferroni_significant_domain_29 %in% TRUE),
  n_primary_eligible = sum(QC_primary_eligible %in% TRUE),
  n_significant_domain_29_primary = sum(Bonferroni_significant_domain_29_primary %in% TRUE)
), by = domain]
summary <- rbind(
  data.table(
    metric = c("n_custom_annotations", "n_raw_domain_rows", "n_finite_raw_domain_contrasts", "n_significant_raw_domain_29", "n_significant_raw_domain_29_primary", "n_warning_rows", "n_error_rows", "raw_alpha_domain_29"),
    value = c(nrow(custom), nrow(raw_table), sum(raw_table$finite_contrast), sum(raw_table$Bonferroni_significant_domain_29), sum(raw_table$Bonferroni_significant_domain_29_primary), sum(raw_table$Annotation_model_warning), sum(raw_table$Annotation_model_error), 0.05 / 29)
  ),
  rbindlist(lapply(seq_len(nrow(domain_summary)), function(i) {
    d <- domain_summary[i]
    data.table(
      metric = c(paste0("n_", d$domain, "_rows"), paste0("n_", d$domain, "_finite"), paste0("n_", d$domain, "_significant_domain_29"), paste0("n_", d$domain, "_primary_eligible"), paste0("n_", d$domain, "_significant_domain_29_primary")),
      value = c(d$n_rows, d$n_finite, d$n_significant_domain_29, d$n_primary_eligible, d$n_significant_domain_29_primary)
    )
  }))
)
fwrite(summary, file.path(out_dir, "D6B_raw_annotation_covariance_summary.tsv"), sep = "\t", na = "NA", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "D6B_raw_covariance_reproduction_sessionInfo.txt"))
cat("RAW_COVARIANCE_DONE annotations=", nrow(custom), " domain_rows=", nrow(raw_table), " finite=", sum(raw_table$finite_contrast), " sig29=", sum(raw_table$Bonferroni_significant_domain_29), " primary_sig29=", sum(raw_table$Bonferroni_significant_domain_29_primary), "\n", sep = "")

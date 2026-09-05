#!/usr/bin/env Rscript

# D5i: primary joint full-sample Genomic SEM downstream analysis.
# PGC and FinnGen SCZ are entered together as external observed traits.
# The three-domain measurement syntax is shared, and uncertainty for the
# standardized target-domain correlations and their differences is propagated
# through the joint LDSC S/V object using GenomicSEM's sandwich construction.

root <- "/public/home/zhengxiaoyang/phenotype_ascertainment_pilot"
base <- file.path(root, "results/d5i_joint_fullsample")
input_path <- file.path(root, "results/d5g_joint_fullsample/joint_covstruc.rds")
dir.create(base, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("missing joint LDSC covariance")

.libPaths(c(file.path(root, "software/Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(lavaan)
  library(GenomicSEM)
})

cov <- readRDS(input_path)
all_traits <- c("SCZ_PGC", "SCZ_FG", "MDD", "BIP", "PTSD", "ANX", "ADHD", "ASD", "OCD", "AN", "AUD", "CUD")
analysis_traits <- c("SCZ_PGC", "SCZ_FG", "MDD", "PTSD", "ANX", "ADHD", "OCD", "AN", "AUD", "CUD")
indicator_traits <- c("MDD", "PTSD", "ANX", "ADHD", "OCD", "AN", "AUD", "CUD")
factors <- c("F_COMP", "F_INT", "F_SUD")
target_names <- c("SCZ_PGC", "SCZ_FG")
if (!identical(colnames(cov$S), all_traits)) stop("unexpected joint S order")
if (!all(dim(cov$S) == c(length(all_traits), length(all_traits)))) stop("unexpected S dimension")
if (!all(dim(cov$V) == c(choose(length(all_traits) + 1L, 2L), choose(length(all_traits) + 1L, 2L)))) stop("unexpected V dimension")
if (!all(is.finite(cov$S)) || !all(is.finite(cov$V))) stop("joint S/V contains non-finite values")

build_model <- function(floor_value = 1e-4) {
  paste(c(
    "F_COMP =~ 1*OCD + l_comp_an*AN + l_comp_anx*ANX",
    "F_INT  =~ 1*MDD + l_int_ptsd*PTSD + l_int_anx*ANX",
    "F_SUD  =~ 1*AUD + l_sud_cud*CUD + l_sud_adhd*ADHD",
    "F_COMP ~~ F_INT + F_SUD",
    "F_INT ~~ F_SUD",
    "F_COMP ~~ v_comp*F_COMP",
    "F_INT  ~~ v_int*F_INT",
    "F_SUD  ~~ v_sud*F_SUD",
    paste0(indicator_traits, " ~~ resid_", indicator_traits, "*", indicator_traits),
    if (!is.null(floor_value)) paste0("resid_OCD > ", format(floor_value, scientific = FALSE, trim = TRUE)),
    "SCZ_PGC ~~ v_pgc*SCZ_PGC + c_pgc_fg*SCZ_FG + c_pgc_comp*F_COMP + c_pgc_int*F_INT + c_pgc_sud*F_SUD",
    "SCZ_FG  ~~ v_fg*SCZ_FG + c_fg_comp*F_COMP + c_fg_int*F_INT + c_fg_sud*F_SUD",
    "r_pgc_comp := c_pgc_comp/sqrt(v_pgc*v_comp)",
    "r_pgc_int  := c_pgc_int/sqrt(v_pgc*v_int)",
    "r_pgc_sud  := c_pgc_sud/sqrt(v_pgc*v_sud)",
    "r_fg_comp  := c_fg_comp/sqrt(v_fg*v_comp)",
    "r_fg_int   := c_fg_int/sqrt(v_fg*v_int)",
    "r_fg_sud   := c_fg_sud/sqrt(v_fg*v_sud)",
    "d_comp := r_fg_comp-r_pgc_comp",
    "d_int  := r_fg_int-r_pgc_int",
    "d_sud  := r_fg_sud-r_pgc_sud"
  ), collapse = "\n")
}

run_usermodel <- function(model_text) {
  warnings_seen <- character()
  printed <- character()
  ans <- NULL
  printed <- capture.output(ans <- tryCatch(withCallingHandlers(
    usermodel(covstruc = list(V = cov$V, S = cov$S), estimation = "DWLS",
      model = model_text, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE,
      fix_resid = FALSE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e))
  list(fit = ans, warnings = unique(warnings_seen), printed = printed)
}

write_run <- function(run, model_id) {
  out <- file.path(base, model_id)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(run$printed, file.path(out, "diagnostic_output.txt"))
  writeLines(run$warnings, file.path(out, "warnings.txt"))
  saveRDS(run$fit, file.path(out, "fit.rds"))
  out
}

finite_all <- function(x) all(is.finite(suppressWarnings(as.numeric(x))))
finite_free <- function(x, p) {
  z <- suppressWarnings(as.numeric(x))
  required <- !(p$op == "=~" & is.na(z))
  any(required) && all(is.finite(z[required]))
}

# Reconstruct the sandwich covariance used by usermodel() so that a joint
# Wald statistic can use the full covariance among the three defined deltas.
defined_sandwich <- function(model_text) {
  V_LD <- as.matrix(cov$V)
  S_LD <- as.matrix(cov$S)
  S_names <- colnames(S_LD)
  y <- expand.grid(S_names, S_names)
  y <- y[!duplicated(apply(y, 1, function(x) paste(sort(x), collapse = ""))), ]
  V_names <- paste(y$Var1, y$Var2, sep = " ")
  colnames(V_LD) <- V_names
  rownames(V_LD) <- V_names

  remove2 <- integer()
  for (i in seq_along(S_names)) {
    if (!grepl(paste0("\\b", S_names[[i]], "\\b"), model_text)) remove2 <- c(remove2, i)
  }
  if (length(remove2)) {
    keep <- setdiff(seq_along(S_names), remove2)
    keep_names <- S_names[keep]
    keep_pairs <- expand.grid(keep_names, keep_names)
    keep_pairs <- keep_pairs[!duplicated(apply(keep_pairs, 1, function(x) paste(sort(x), collapse = ""))), ]
    keep_vnames <- paste(keep_pairs$Var1, keep_pairs$Var2, sep = " ")
    V_LD <- V_LD[keep_vnames, keep_vnames, drop = FALSE]
    S_LD <- S_LD[keep, keep, drop = FALSE]
  }
  k <- ncol(S_LD)
  z <- choose(k + 1L, 2L)
  eig_s <- eigen(S_LD, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig_s) <= 0) S_LD <- as.matrix(Matrix::nearPD(S_LD, corr = FALSE)$mat)
  eig_v <- eigen(V_LD, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig_v) <= 0) V_LD <- as.matrix(Matrix::nearPD(V_LD, corr = FALSE)$mat)

  reorder_fit <- lavaan::sem(model_text, sample.cov = S_LD, estimator = "DWLS",
    WLS.V = solve(V_LD), se = "standard", sample.nobs = 2,
    std.lv = FALSE, optim.dx.tol = 0.01, optim.force.converged = TRUE,
    control = list(iter.max = 1), warn = FALSE)
  order <- GenomicSEM:::.rearrange(k = k, fit = reorder_fit, names = rownames(S_LD))
  V_reorder <- V_LD[order, order, drop = FALSE]
  W_reorder <- diag(z)
  diag(W_reorder) <- diag(V_reorder)
  W_reorder <- solve(W_reorder)

  primary_fit <- lavaan::sem(model_text, sample.cov = S_LD, estimator = "DWLS",
    WLS.V = W_reorder, se = "standard", sample.nobs = 2,
    std.lv = FALSE, optim.dx.tol = 0.01, optim.force.converged = TRUE,
    warn = FALSE)
  delta <- lavInspect(primary_fit, "delta")
  W <- lavInspect(primary_fit, "WLS.V")
  bread <- solve(t(delta) %*% W %*% delta)
  lettuce <- W %*% delta
  Ohtt <- bread %*% t(lettuce) %*% V_reorder %*% lettuce %*% bread

  pt <- parTable(primary_fit)
  def <- pt[pt$op == ":=", , drop = FALSE]
  if (!nrow(def)) stop("no defined parameters in lavaan model")
  func <- primary_fit@Model@def.function
  x <- lavaan:::lav_model_get_parameters(primary_fit@Model, type = "free")
  Jac <- lavaan:::lav_func_jacobian_complex(func = func, x = x)
  col_use <- which(colSums(abs(Jac)) > 0)
  Jac_sub <- Jac[, col_use, drop = FALSE]
  V_sub <- Ohtt[col_use, col_use, drop = FALSE]
  V_sub <- 0.5 * (V_sub + t(V_sub))
  eig <- eigen(V_sub, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || min(eig) < -1e-10) {
    ee <- eigen(V_sub, symmetric = TRUE)
    V_sub <- ee$vectors %*% (pmax(ee$values, 0) * t(ee$vectors))
    V_sub <- 0.5 * (V_sub + t(V_sub))
  }
  V_def <- Jac_sub %*% V_sub %*% t(Jac_sub)
  dimnames(V_def) <- list(def$lhs, def$lhs)
  optim_info <- tryCatch(lavInspect(primary_fit, "optim"), error = function(e) list())
  list(primary_fit = primary_fit, V_reorder = V_reorder, Ohtt = Ohtt,
    V_def = V_def, defined_table = def, bread = bread, optim = optim_info)
}

optimizer_audit <- function(sandwich, model_id, floor_value = NULL) {
  if (inherits(sandwich, "error") || is.null(sandwich)) {
    return(data.table(model_id = model_id, lavaan_optim_converged = FALSE,
      max_abs_dx_all = NA_real_, max_abs_dx_excluding_active_boundary_parameters = NA_real_,
      parameter_with_max_dx = "", OCD_residual_estimate = NA_real_,
      OCD_constraint_active = NA, error = if (inherits(sandwich, "error")) conditionMessage(sandwich) else "sandwich reconstruction unavailable"))
  }
  pt <- parTable(sandwich$primary_fit)
  free <- pt[pt$free > 0 & !duplicated(pt$free), , drop = FALSE]
  free <- free[order(free$free), , drop = FALSE]
  dx <- as.numeric(sandwich$optim$dx)
  max_idx <- if (length(dx)) which.max(abs(dx)) else NA_integer_
  label_for <- function(i) {
    if (!is.finite(i) || i < 1L || i > nrow(free)) return("")
    z <- free[i, , drop = FALSE]
    paste(z$lhs, z$op, z$rhs, if (nzchar(z$label)) paste0("[", z$label, "]") else "")
  }
  ocd <- free[free$label == "resid_OCD", , drop = FALSE]
  ocd_est <- if (nrow(ocd)) as.numeric(ocd$est[[1L]]) else NA_real_
  active <- if (nrow(ocd) && !is.null(floor_value)) is.finite(ocd_est) && ocd_est <= floor_value * 1.01 else FALSE
  active_idx <- if (active) as.integer(ocd$free[[1L]]) else integer()
  non_active <- setdiff(seq_along(dx), active_idx)
  non_active_max <- if (length(non_active)) max(abs(dx[non_active])) else NA_real_
  data.table(model_id = model_id,
    lavaan_optim_converged = isTRUE(sandwich$optim$converged),
    max_abs_dx_all = if (length(dx)) max(abs(dx)) else NA_real_,
    max_abs_dx_excluding_active_boundary_parameters = non_active_max,
    parameter_with_max_dx = label_for(max_idx), OCD_residual_estimate = ocd_est,
    OCD_constraint_active = active, error = "")
}

assess <- function(run, model_id, floor_value) {
  out <- write_run(run, model_id)
  fit <- run$fit
  if (inherits(fit, "error") || is.null(fit) || is.null(fit$results) || !nrow(as.data.frame(fit$results))) {
    row <- data.table(model_id = model_id, boundary_floor = floor_value,
      returned_parameter_table = FALSE, gate_pass = FALSE,
      error = if (inherits(fit, "error")) conditionMessage(fit) else "No interpretable parameter table")
    fwrite(row, file.path(out, "gate.tsv"), sep = "\t")
    return(list(row = row, sandwich = NULL))
  }

  p <- as.data.table(fit$results)
  fwrite(p, file.path(out, "parameters.tsv"), sep = "\t")
  fwrite(as.data.table(fit$modelfit), file.path(out, "fit_metrics.tsv"), sep = "\t")
  est_col <- intersect(c("Unstand_Est", "est"), names(p))[[1L]]
  se_col <- intersect(c("Unstand_SE", "se"), names(p))[[1L]]
  std_col <- intersect(c("STD_All", "STD_Genotype"), names(p))[[1L]]
  defs <- p[op == ":=" & lhs %chin% c("r_pgc_comp", "r_pgc_int", "r_pgc_sud",
    "r_fg_comp", "r_fg_int", "r_fg_sud", "d_comp", "d_int", "d_sud")]
  defs[, `:=`(estimate = suppressWarnings(as.numeric(get(est_col))),
    SE = suppressWarnings(as.numeric(get(se_col))))]
  fwrite(defs, file.path(out, "defined_parameters.tsv"), sep = "\t")
  residuals <- p[op == "~~" & lhs == rhs & lhs %chin% indicator_traits]
  latent_var <- p[op == "~~" & lhs == rhs & lhs %chin% factors]
  neg <- residuals[is.finite(suppressWarnings(as.numeric(get(std_col)))) &
    suppressWarnings(as.numeric(get(std_col))) < 0, lhs]
  mf <- as.data.table(fit$modelfit)
  cfi <- as.numeric(mf$CFI[[1L]])
  srmr <- as.numeric(mf$SRMR[[1L]])
  sandwich <- tryCatch(defined_sandwich(build_model(floor_value)), error = function(e) e)
  sandwich_ok <- !inherits(sandwich, "error") && all(is.finite(sandwich$V_def))
  if (sandwich_ok) saveRDS(sandwich, file.path(out, "sandwich_reconstruction.rds"))
  d_names <- c("d_comp", "d_int", "d_sud")
  d_idx <- if (sandwich_ok) match(d_names, rownames(sandwich$V_def)) else integer()
  d_est <- defs[match(d_names, lhs), estimate]
  d_se_rebuilt <- if (sandwich_ok) sqrt(pmax(diag(sandwich$V_def)[d_idx], 0)) else rep(NA_real_, 3L)
  d_cov <- if (sandwich_ok) sandwich$V_def[d_idx, d_idx, drop = FALSE] else matrix(NA_real_, 3L, 3L)
  dimnames(d_cov) <- list(d_names, d_names)
  rank_d <- if (sandwich_ok) qr(d_cov, tol = 1e-10)$rank else NA_integer_
  omnibus_w <- NA_real_; omnibus_p <- NA_real_
  if (sandwich_ok && rank_d == 3L && all(is.finite(d_est)) && all(is.finite(d_cov))) {
    omnibus_w <- as.numeric(t(d_est) %*% solve(d_cov, d_est))
    omnibus_p <- pchisq(omnibus_w, df = 3L, lower.tail = FALSE)
  }
  if (sandwich_ok) {
    fwrite(as.data.table(sandwich$V_def, keep.rownames = "parameter"),
      file.path(out, "defined_parameter_covariance.tsv"), sep = "\t")
    audit <- data.table(
      sandwich_finite = all(is.finite(sandwich$Ohtt)),
      defined_covariance_finite = all(is.finite(sandwich$V_def)),
      defined_covariance_rank = rank_d,
      lavaan_optim_converged = isTRUE(sandwich$optim$converged),
      lavaan_max_abs_dx = if (length(sandwich$optim$dx)) max(abs(sandwich$optim$dx)) else NA_real_,
      d_cov_min_eigen = min(eigen((d_cov + t(d_cov))/2, symmetric = TRUE, only.values = TRUE)$values),
      d_omnibus_W = omnibus_w, d_omnibus_df = rank_d, d_omnibus_P = omnibus_p,
      d_comp_SE_rebuilt = d_se_rebuilt[[1L]], d_int_SE_rebuilt = d_se_rebuilt[[2L]],
      d_sud_SE_rebuilt = d_se_rebuilt[[3L]]
    )
    fwrite(audit, file.path(out, "sandwich_audit.tsv"), sep = "\t")
  }
  hard <- finite_all(p[[est_col]]) && finite_all(p[[std_col]]) && finite_free(p[[se_col]], p) &&
    nrow(residuals) == length(indicator_traits) &&
    all(suppressWarnings(as.numeric(residuals[[est_col]])) >= -1e-10) &&
    nrow(latent_var) == length(factors) &&
    all(suppressWarnings(as.numeric(latent_var[[est_col]])) > 0) &&
    nrow(defs) == 9L && finite_all(defs$estimate) && finite_all(defs$SE) && sandwich_ok &&
    isTRUE(sandwich$optim$converged)
  fit_ok <- is.finite(cfi) && cfi >= 0.90 && is.finite(srmr) && srmr <= 0.10
  row <- data.table(model_id = model_id, boundary_floor = floor_value,
    returned_parameter_table = TRUE, CFI = cfi, SRMR = srmr,
    negative_residuals = paste(neg, collapse = ","),
    hard_admissibility = hard, acceptable_fit = fit_ok,
    defined_covariance_rank = rank_d, omnibus_W = omnibus_w,
    omnibus_df = if (is.finite(omnibus_w)) rank_d else NA_integer_,
    omnibus_P = omnibus_p, gate_pass = hard && fit_ok, error = "")
  fwrite(row, file.path(out, "gate.tsv"), sep = "\t")
  list(row = row, sandwich = sandwich, defs = defs)
}

primary <- assess(run_usermodel(build_model(1e-4)), "primary_boundary_0.0001", 1e-4)
sensitivity <- assess(run_usermodel(build_model(1e-3)), "sensitivity_boundary_0.001", 1e-3)
summary <- rbindlist(list(primary$row, sensitivity$row), fill = TRUE)
fwrite(summary, file.path(base, "D5i_fit_summary.tsv"), sep = "\t")

# Convergence audit: retain the unconstrained diagnostic, while using the
# official GenomicSEM/lavaan convergence state for the constrained model. The
# raw derivative at an active inequality boundary is reported, not gated.
free_model <- build_model(NULL)
free_run <- run_usermodel(free_model)
free_out <- write_run(free_run, "unconstrained_diagnostic")
free_sandwich <- tryCatch(defined_sandwich(free_model), error = function(e) e)
constrained_sandwich <- primary$sandwich
conv_audit <- rbindlist(list(
  optimizer_audit(free_sandwich, "unconstrained", NULL),
  optimizer_audit(constrained_sandwich, "constrained_boundary_0.0001", 1e-4)
), fill = TRUE)
if (inherits(free_run$fit, "error") || is.null(free_run$fit) || is.null(free_run$fit$results)) {
  conv_audit[, `:=`(standardized_model_available = FALSE, sandwich_vcov_available = !inherits(free_sandwich, "error"))]
} else {
  fp <- as.data.table(free_run$fit$results)
  conv_audit[model_id == "unconstrained", `:=`(
    standardized_model_available = "STD_All" %in% names(fp),
    sandwich_vcov_available = !inherits(free_sandwich, "error") && all(is.finite(free_sandwich$Ohtt))
  )]
}
cp <- primary$row
conv_audit[model_id == "constrained_boundary_0.0001", `:=`(
  standardized_model_available = isTRUE(primary$row$returned_parameter_table[[1L]]),
  sandwich_vcov_available = !is.null(primary$sandwich) && all(is.finite(primary$sandwich$Ohtt)),
  latent_variances_admissible = isTRUE(primary$row$hard_admissibility[[1L]]),
  latent_correlations_admissible = isTRUE(primary$row$hard_admissibility[[1L]])
)]
fwrite(conv_audit, file.path(base, "convergence_audit.tsv"), sep = "\t")

# Direct observed-trait comparator differences use the same joint LDSC V.
v_names <- paste((function() {
  y <- expand.grid(all_traits, all_traits)
  y[!duplicated(apply(y, 1, function(x) paste(sort(x), collapse = ""))), ]
})()$Var1, (function() {
  y <- expand.grid(all_traits, all_traits)
  y[!duplicated(apply(y, 1, function(x) paste(sort(x), collapse = ""))), ]
})()$Var2, sep = " ")
pair_index <- function(a, b) {
  hit <- match(c(paste(a, b), paste(b, a)), v_names)
  hit[which(!is.na(hit))[1L]]
}
direct_rg_grad <- function(a, b, S) {
  saa <- S[a, a]; sbb <- S[b, b]; sab <- S[a, b]
  rg <- sab / sqrt(saa * sbb)
  g <- numeric(length(v_names))
  g[pair_index(a, b)] <- 1 / sqrt(saa * sbb)
  g[pair_index(a, a)] <- -0.5 * rg / saa
  g[pair_index(b, b)] <- -0.5 * rg / sbb
  list(rg = as.numeric(rg), grad = g)
}
comparator_rows <- list()
for (trait in c("BIP", "ASD")) {
  pg <- direct_rg_grad("SCZ_PGC", trait, cov$S)
  fg <- direct_rg_grad("SCZ_FG", trait, cov$S)
  grad <- fg$grad - pg$grad
  vv <- as.numeric(t(grad) %*% cov$V %*% grad)
  se <- if (is.finite(vv) && vv >= 0) sqrt(vv) else NA_real_
  delta <- fg$rg - pg$rg
  comparator_rows[[length(comparator_rows) + 1L]] <- data.table(
    comparator = trait, rg_pgc = pg$rg, rg_finngen = fg$rg, delta = delta,
    SE = se, Z = delta / se, P = 2 * pnorm(abs(delta / se), lower.tail = FALSE),
    variance_from_joint_V = vv
  )
}
fwrite(rbindlist(comparator_rows), file.path(base, "observed_comparator_deltas.tsv"), sep = "\t")

primary_ok <- isTRUE(primary$row$gate_pass[[1L]])
sensitivity_ok <- isTRUE(sensitivity$row$gate_pass[[1L]])
status <- c(
  "analysis=D5i exploratory joint full-sample Genomic SEM transport analysis",
  "measurement_structure=F_COMP,F_INT,F_SUD; SCZ_PGC and SCZ_FG are external observed traits",
  "joint_input=existing D5g 12-trait full-sample LDSC S/V; D5i primary model uses both SCZ targets and 8 indicators",
  "inference_source=joint LDSC S/V with GenomicSEM sandwich covariance and defined-parameter delta method",
  "primary_endpoint=three-dimensional standardized target-domain correlation profile difference",
  "secondary=BIP and ASD direct observed-trait comparator differences from the same joint S/V",
  "additional_block_refit=failed D5h stability sensitivity retained as diagnostic and not used as a gate",
  "boundary_caveat=OCD residual reached the prespecified lower boundary; sandwich reconstruction recorded an active-boundary optimizer warning",
  "interpretation_scope=exploratory full-sample boundary-constrained Wald/delta inference; report the boundary estimate explicitly",
  paste0("primary_gate=", ifelse(primary_ok, "PASS", "STOP")),
  paste0("boundary_sensitivity_gate=", ifelse(sensitivity_ok, "PASS", "STOP")),
  paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
)
if (!primary_ok) {
  status <- c(status, "formal_inference=STOP", "stop_reason=primary joint full-sample model failed admissibility or fit gate")
} else {
  status <- c(status, "formal_inference=PASS", "formal_inference_scope=full-sample joint LDSC S/V sandwich inference; no SEM block-refit gate")
}
writeLines(status, file.path(base, "D5i_STATUS.txt"))
cat("D5I_JOINT_MODEL_COMPLETE\n")

#!/usr/bin/env Rscript

# D6A: official BaselineLD v2.2 Stratified Genomic SEM analysis.
# The model syntax follows D5i. PGC and FinnGen SCZ are jointly entered as
# observed target traits, and the three D5i latent domains are retained.
# Uncertainty is obtained from s_ldsc()'s 200-block LDSC jackknife and the
# sandwich covariance used by enrich(); no SEM block-refit gate is applied.

root <- "/public/home/zhengxiaoyang/phenotype_ascertainment_pilot"
out_dir <- file.path(root, "results/d6a_baseline_stratified_gsem_20260905")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "analysis=D6A official BaselineLD v2.2 Stratified Genomic SEM",
  "status=RUNNING",
  paste0("started_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(out_dir, "D6A_STATUS.txt"))

options(timeout = 1800)
.libPaths(c(file.path(root, "software/Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

traits <- c(
  SCZ_PGC = file.path(root, "results/p11_scz_neff_audit/input/SCZ_PGC2022_Nx2.sumstats.gz"),
  SCZ_FG  = file.path(root, "results/p1/ldsc_input/SCZ_FINNGEN_R13.sumstats.gz"),
  MDD     = file.path(root, "results/p1/ldsc_input/MDD_CLIN_PGC2025.sumstats.gz"),
  PTSD    = "/public/home/zhengxiaoyang/BPD_transdiagnostic_genetics/analysis/ldsc/munged/PTSD.sumstats.gz",
  ANX     = "/public/home/zhengxiaoyang/BPD_transdiagnostic_genetics/analysis/ldsc/munged/ANX.sumstats.gz",
  ADHD    = "/public/home/zhengxiaoyang/BPD_transdiagnostic_genetics/analysis/ldsc/munged/ADHD.sumstats.gz",
  OCD     = file.path(root, "results/p25/aux_munged/OCD_2025.sumstats.gz"),
  AN      = file.path(root, "results/p3/aux_munged/AN.sumstats.gz"),
  AUD     = "/public/home/zhengxiaoyang/BPD_transdiagnostic_genetics/analysis/ldsc/munged/AUD.sumstats.gz",
  CUD     = file.path(root, "results/p25/aux_munged/CUD_2023_EUR.sumstats.gz")
)

baseline_ld <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/09_single_cell/dice_ldsc_seg/reference/"
wld <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/13_ldsc_seg_official/reference/1000G_Phase3_weights_hm3_no_MHC/"
frq <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/09_single_cell/dice_ldsc_seg/reference/1000G_Phase3_frq/"
reuse_existing_s_ldsc <- identical(Sys.getenv("D6A_REUSE_LDSC"), "1")

header_qc <- rbindlist(lapply(names(traits), function(nm) {
  path <- unname(traits[[nm]])
  hdr <- fread(cmd = paste("zcat", shQuote(path)), nrows = 0L, showProgress = FALSE)
  data.table(
    trait = nm, path = path, exists = file.exists(path),
    has_SNP = "SNP" %in% names(hdr), has_A1 = "A1" %in% names(hdr),
    has_A2 = "A2" %in% names(hdr), has_Z = "Z" %in% names(hdr),
    has_N = "N" %in% names(hdr), columns = paste(names(hdr), collapse = ",")
  )
}), fill = TRUE)
fwrite(header_qc, file.path(out_dir, "input_header_qc.tsv"), sep = "\t")
if (!all(header_qc$exists & header_qc$has_SNP & header_qc$has_A1 & header_qc$has_A2 & header_qc$has_Z & header_qc$has_N)) {
  stop("One or more D5i summary-statistic inputs failed existence/header QC")
}

ref_qc <- data.table(
  reference = c("baselineLD_v2.2", "weights_hm3_no_MHC", "frq"),
  path = c(baseline_ld, wld, frq),
  exists = c(dir.exists(baseline_ld), dir.exists(wld), dir.exists(frq))
)
ref_qc[, n_annot := c(
  sum(grepl("baselineLD\\.[0-9]+\\.annot\\.gz$", list.files(baseline_ld))),
  NA_integer_, NA_integer_
)]
ref_qc[, n_ldscore := c(
  sum(grepl("baselineLD\\.[0-9]+\\.l2\\.ldscore\\.gz$", list.files(baseline_ld))),
  sum(grepl("weights\\.hm3_noMHC\\.[0-9]+\\.l2\\.ldscore\\.gz$", list.files(wld))),
  NA_integer_
)]
ref_qc[, n_frq := c(NA_integer_, NA_integer_, sum(grepl("1000G\\.EUR\\.QC\\.[0-9]+\\.frq$", list.files(frq))))]
fwrite(ref_qc, file.path(out_dir, "reference_qc.tsv"), sep = "\t")
if (!all(ref_qc$exists) || ref_qc[reference == "baselineLD_v2.2", n_annot] != 22L ||
    ref_qc[reference == "baselineLD_v2.2", n_ldscore] != 22L ||
    ref_qc[reference == "weights_hm3_no_MHC", n_ldscore] != 22L ||
    ref_qc[reference == "frq", n_frq] != 22L) {
  stop("BaselineLD v2.2, weights, or frequency reference files are incomplete")
}

model <- paste(c(
  "F_COMP =~ 1*OCD + l_comp_an*AN + l_comp_anx*ANX",
  "F_INT  =~ 1*MDD + l_int_ptsd*PTSD + l_int_anx*ANX",
  "F_SUD  =~ 1*AUD + l_sud_cud*CUD + l_sud_adhd*ADHD",
  "F_COMP ~~ F_INT + F_SUD",
  "F_INT ~~ F_SUD",
  "F_COMP ~~ v_comp*F_COMP",
  "F_INT  ~~ v_int*F_INT",
  "F_SUD  ~~ v_sud*F_SUD",
  paste0(c("MDD", "PTSD", "ANX", "ADHD", "OCD", "AN", "AUD", "CUD"),
         " ~~ resid_", c("MDD", "PTSD", "ANX", "ADHD", "OCD", "AN", "AUD", "CUD"),
         "*", c("MDD", "PTSD", "ANX", "ADHD", "OCD", "AN", "AUD", "CUD")),
  "resid_OCD > 0.0001",
  "SCZ_PGC ~~ v_pgc*SCZ_PGC + c_pgc_fg*SCZ_FG + c_pgc_comp*F_COMP + c_pgc_int*F_INT + c_pgc_sud*F_SUD",
  "SCZ_FG  ~~ v_fg*SCZ_FG + c_fg_comp*F_COMP + c_fg_int*F_INT + c_fg_sud*F_SUD"
), collapse = "\n")

params <- c(
  # lavaan canonicalizes observed--latent covariance rows with the latent
  # variable on the left-hand side; enrich() matches this exact ordering.
  "F_COMP~~SCZ_PGC", "F_INT~~SCZ_PGC", "F_SUD~~SCZ_PGC",
  "F_COMP~~SCZ_FG",  "F_INT~~SCZ_FG",  "F_SUD~~SCZ_FG"
)

config <- c(
  paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("GenomicSEM_version=", as.character(packageVersion("GenomicSEM"))),
  paste0("R_version=", R.version.string),
  "analysis=D6A official BaselineLD v2.2 Stratified Genomic SEM",
  "model_source=D5i three-domain joint external-trait model",
  paste0("traits=", paste(names(traits), collapse = ",")),
  paste0("parameters=", paste(params, collapse = ",")),
  "sample_prev=NA for all traits; effective-N sumstats are used",
  "n.blocks=200",
  "tau=FALSE",
  "fix=regressions",
  "std.lv=FALSE",
  "base=TRUE",
  "rm_flank=TRUE",
  "boundary=resid_OCD > 0.0001; inherited from D5i; no new boundary",
  paste0("baseline_ld=", baseline_ld),
  paste0("wld=", wld),
  paste0("frq=", frq),
  "multiple_testing=Bonferroni across six target parameters and analyzed non-base BaselineLD annotations",
  "inference_source=LDSC 200-block jackknife V and GenomicSEM enrich sandwich covariance",
  "additional_SE_gate=none",
  paste0("reuse_existing_s_ldsc=", reuse_existing_s_ldsc)
)
writeLines(config, file.path(out_dir, "D6A_CONFIG.txt"))
writeLines(model, file.path(out_dir, "D6A_model_syntax.txt"))
fwrite(data.table(trait = names(traits), path = unname(traits)), file.path(out_dir, "trait_manifest.tsv"), sep = "\t")

fail <- function(e) {
  writeLines(c(config, "status=ERROR", paste0("error=", conditionMessage(e)),
               paste0("finished_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))),
             file.path(out_dir, "D6A_STATUS.txt"))
  stop(e)
}

tryCatch({
  if (reuse_existing_s_ldsc) {
    cached_s <- file.path(out_dir, "D6A_S_LDSCOutput.rds")
    if (!file.exists(cached_s)) stop("D6A_REUSE_LDSC=1 but cached S/LDSC object is missing")
    cat("D6A_S_LDSC_REUSED\n", flush = TRUE)
    S_LDSCOutput <- readRDS(cached_s)
  } else {
    cat("D6A_S_LDSC_START\n", flush = TRUE)
    S_LDSCOutput <- s_ldsc(
      traits = unname(traits),
      sample.prev = rep(NA_real_, length(traits)),
      population.prev = rep(NA_real_, length(traits)),
      ld = baseline_ld,
      wld = wld,
      frq = frq,
      trait.names = names(traits),
      n.blocks = 200L,
      ldsc.log = file.path(out_dir, "D6A_s_ldsc"),
      exclude_cont = TRUE
    )
    saveRDS(S_LDSCOutput, file.path(out_dir, "D6A_S_LDSCOutput.rds"), compress = FALSE)
  }
  if (is.null(S_LDSCOutput$S) || is.null(S_LDSCOutput$V)) stop("s_ldsc returned no S/V")
  cat("D6A_S_LDSC_DONE\n", flush = TRUE)

  n_s <- length(S_LDSCOutput$S)
  if (n_s < 2L) stop("s_ldsc returned fewer than two annotation matrices")
  cat("D6A_ENRICH_START\n", flush = TRUE)
  E <- enrich(
    s_covstruc = S_LDSCOutput,
    model = model,
    params = params,
    fix = "regressions",
    std.lv = FALSE,
    rm_flank = TRUE,
    tau = FALSE,
    base = TRUE
  )
  saveRDS(E, file.path(out_dir, "D6A_enrich_raw.rds"), compress = FALSE)
  cat("D6A_ENRICH_DONE\n", flush = TRUE)

  pieces <- lapply(seq_along(E), function(i) {
    z <- as.data.table(E[[i]])
    z[, target_parameter := if (i <= length(params)) params[[i]] else "BASELINE"]
    z
  })
  enrichment <- rbindlist(pieces, fill = TRUE)
  n_base <- sum(enrichment$target_parameter == "BASELINE")
  tested <- enrichment[target_parameter != "BASELINE" & !is.na(Annotation) &
    !tolower(as.character(Annotation)) %in% c("base", "baseline")]
  n_annotation_tests <- length(unique(tested$Annotation))
  if (n_annotation_tests < 1L) stop("No non-base annotation rows were returned by enrich")
  alpha <- 0.05 / (n_annotation_tests * length(params))
  enrichment[, Bonferroni_alpha := alpha]
  enrichment[, Bonferroni_significant := !is.na(Enrichment_p_value) & Enrichment_p_value < alpha]
  enrichment[, Annotation_model_warning := !is.na(Warning) & trimws(as.character(Warning)) != "0"]
  enrichment[, Annotation_model_error := !is.na(Error) & trimws(as.character(Error)) != "0"]
  enrichment[, Bonferroni_significant_no_warning := Bonferroni_significant &
    !Annotation_model_warning & !Annotation_model_error]
  enrichment[, n_annotations_in_family := n_annotation_tests]
  enrichment[, n_parameters_in_family := length(params)]
  fwrite(enrichment, file.path(out_dir, "D6A_enrichment_long.tsv"), sep = "\t", na = "NA")

  sig <- enrichment[target_parameter != "BASELINE" & Bonferroni_significant == TRUE]
  sig_no_warning <- enrichment[target_parameter != "BASELINE" & Bonferroni_significant_no_warning == TRUE]
  fwrite(sig, file.path(out_dir, "D6A_enrichment_bonferroni_significant.tsv"), sep = "\t", na = "NA")
  fwrite(sig_no_warning, file.path(out_dir, "D6A_enrichment_bonferroni_significant_no_warning.tsv"), sep = "\t", na = "NA")
  counts <- enrichment[target_parameter != "BASELINE", .(
    n_rows = .N,
    n_finite_enrichment = sum(is.finite(Enrichment)),
    n_finite_p = sum(is.finite(Enrichment_p_value)),
    n_bonferroni_significant = sum(Bonferroni_significant, na.rm = TRUE),
    n_bonferroni_significant_no_warning = sum(Bonferroni_significant_no_warning, na.rm = TRUE),
    n_annotation_warnings = sum(Annotation_model_warning, na.rm = TRUE),
    n_annotation_errors = sum(Annotation_model_error, na.rm = TRUE),
    min_p = suppressWarnings(min(Enrichment_p_value, na.rm = TRUE))
  ), by = target_parameter]
  fwrite(counts, file.path(out_dir, "D6A_parameter_summary.tsv"), sep = "\t", na = "NA")
  fwrite(data.table(
    n_s_ldsc_matrices = n_s,
    n_enrich_list_elements = length(E),
    n_base_rows = n_base,
    n_annotation_tests = n_annotation_tests,
    n_parameters = length(params),
    bonferroni_alpha = alpha,
    n_significant = nrow(sig),
    n_significant_no_warning = nrow(sig_no_warning)
  ), file.path(out_dir, "D6A_family_summary.tsv"), sep = "\t", na = "NA")
  writeLines(c(config, paste0("status=COMPLETED"), paste0("n_s_ldsc_matrices=", n_s),
               paste0("n_enrich_list_elements=", length(E)),
               paste0("n_annotation_tests=", n_annotation_tests),
               paste0("bonferroni_alpha=", format(alpha, scientific = TRUE)),
               paste0("n_bonferroni_significant=", nrow(sig)),
               paste0("n_bonferroni_significant_no_warning=", nrow(sig_no_warning)),
               paste0("finished_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))),
             file.path(out_dir, "D6A_STATUS.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  cat("D6A_COMPLETE\n", flush = TRUE)
}, error = fail)

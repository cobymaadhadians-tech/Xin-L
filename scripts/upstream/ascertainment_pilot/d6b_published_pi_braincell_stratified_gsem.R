#!/usr/bin/env Rscript

# D6B: published PI x brain-cell Stratified Genomic SEM protocol.
#
# The annotation universe follows Supplementary Table 10 of the published
# protocol: 97 BaselineLD v2.2 annotations, 24 selected expression
# annotations, 60 selected Roadmap chromatin annotations, and 29 custom
# pLI-derived PI/cell annotations. The 168 binary annotations define the
# published enrichment-testing family; continuous and flanking BaselineLD
# annotations remain in the conditioning matrix and are removed only by the
# official enrich() settings where appropriate.
#
# Inference uses the official LDSC block-jackknife S/V and enrich() sandwich
# covariance. No second-layer SEM block-refit gate is applied.

root <- "/public/home/zhengxiaoyang/phenotype_ascertainment_pilot"
out_dir <- file.path(root, "results/d6b_published_pi_braincell_stratified_gsem_20260905")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

status_file <- file.path(out_dir, "D6B_STATUS.txt")
started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
writeLines(c(
  "analysis=D6B published PI x brain-cell Stratified Genomic SEM",
  "status=RUNNING",
  paste0("started_utc=", started_utc)
), status_file)

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

baseline_dir <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/09_single_cell/dice_ldsc_seg/reference/"
wld <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/13_ldsc_seg_official/reference/1000G_Phase3_weights_hm3_no_MHC/"
frq <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/09_single_cell/dice_ldsc_seg/reference/1000G_Phase3_frq/"
expr_dir <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/13_ldsc_seg_official/reference/Multi_tissue_gene_expr_1000Gv3_ldscores"
chrom_dir <- "/public/home/zhengxiaoyang/MDD_SIN_analysis/13_ldsc_seg_official/reference/Multi_tissue_chromatin_1000Gv3_ldscores"
custom_dir <- file.path(root, ".codex_work/d6_stratified_gsem_20260905/d6b_annotations/ldscore")
custom_prefix <- file.path(custom_dir, "d6b.")

expression_specs <- data.table(
  annotation_class = "expression",
  published_label = c(
    "Brain Stem", "Whole Brain", "Epidermis", "Fibroblasts",
    "Epithelial Cells", "Embryoid Bodies", "Adipose Subcutaneous",
    "Amygdala", "Cingulate Cortex", "Caudate", "Cerebellar Hemispher",
    "Cerebellum", "Cortex", "Frontal Cortex", "Hippocampus",
    "Hypothalamus", "Nucleus Accumbens", "Putamen", "Spinal Cord",
    "Substantia Nigra", "Endocrine Cells", "Ovary", "Pituitary", "Testis"
  ),
  source_label = c(
    "A08.186.211.132.Brain.Stem", "A08.186.211.Brain",
    "A10.272.497.Epidermis", "A11.329.228.Fibroblasts",
    "A11.436.Epithelial.Cells", "A11.872.190.260.Embryoid.Bodies",
    "Adipose_Subcutaneous", "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_(BA24)",
    "Brain_Caudate_(basal_ganglia)", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_(BA9)",
    "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_(basal_ganglia)",
    "Brain_Putamen_(basal_ganglia)", "Brain_Spinal_cord_(cervical_c-1)",
    "Brain_Substantia_nigra", "A11.382.Endocrine.Cells", "Ovary",
    "Pituitary", "Testis"
  ),
  source_prefix = c(
    "Franke.13.", "Franke.14.", "Franke.31.", "Franke.39.",
    "Franke.30.", "Franke.147.", "GTEx.1.", "GTEx.8.", "GTEx.9.",
    "GTEx.10.", "GTEx.11.", "GTEx.12.", "GTEx.13.", "GTEx.14.",
    "GTEx.15.", "GTEx.16.", "GTEx.17.", "GTEx.18.", "GTEx.19.",
    "GTEx.20.", "Franke.152.", "GTEx.40.", "GTEx.42.", "GTEx.49."
  )
)
expression_specs[, prefix := file.path(expr_dir, source_prefix)]

chromatin_specs <- data.table(
  annotation_class = "chromatin",
  published_label = c(
    "Adipose Nuclei H3K4me1",
    "Angular Gyrus H3K27ac", "Angular Gyrus H3K36me3",
    "Angular Gyrus H3K4me1", "Angular Gyrus H3K4me3",
    "Angular Gyrus H3K9ac",
    "Anterior Caudate H3K27ac", "Anterior Caudate H3K36me3",
    "Anterior Caudate H3K4me1", "Anterior Caudate H3K4me3",
    "Anterior Caudate H3K9ac",
    "Cingulate Gyrus H3K27ac", "Cingulate Gyrus H3K36me3",
    "Cingulate Gyrus H3K4me1", "Cingulate Gyrus H3K4me3",
    "Cingular Gyrus H3K9ac",
    "dLPFC H3K27ac", "dLPFC H3K36me3", "dLPFC H3K4me1",
    "dLPFC H3K4me3", "dLPFC H3K9ac",
    "Germinal Matrix H3K6me3", "Germinal Matrix H3K4me1",
    "Germinal Matrix H3K4me3",
    "Middle Hippocampus H3K27ac", "Middle Hippocampus H3K36me3",
    "Middle Hippocampus H3K4me1", "Middle Hippocampus H3K4me3",
    "Inferior Temporal Lobe H3K27ac", "Interior Temporal Lobe H3K36me3",
    "Inferior Temporal Lobe H3K4me1", "Interior Temporal Lobe H3K4me3",
    "Inferior Temporal Lobe H3K9ac",
    "Substantia Nigra H3K27ac", "Substantia Nigra H3K36me3",
    "Substantia Nigra H3K4me1", "Substantia Nigra H3K4me3",
    "Subtantia Nigra H3K9ac",
    "Fetal Adrenal Gland DNase", "Fetal Adrenal Gland H3K27ac",
    "Fetal Adrenal Gland H3K36me3", "Fetal Adrenal Gland H3K4me1",
    "Fetal Adrenal Gland H3K4me3",
    "Fetal Female Brain Dnase", "Fetal Female Brain H3K36me3",
    "Fetal Female Brain H3K4me1", "Fetal Female Brain H3K4me3",
    "Fetal Male Brain DNAse", "Fetal Male Brain H3K36me3",
    "Fetal Male Brain H3K4me1", "Fetal Male Brain H3K4me3",
    "Foreskin H3K27ac", "Mammry Epithelial Cells H3K27ac",
    "Ovary DNAse", "Ovary H3K27ac", "Ovary H3K36me3",
    "Ovary H3K4me1", "Ovary H3K4me3",
    "Monoculear Peripheral Blood Cells H3K4me1",
    "Rectal Mucosa H3K4me1"
  ),
  source_label = c(
    "Adipose_Nuclei__H3K4me1",
    "Brain_Angular_Gyrus__H3K27ac", "Brain_Angular_Gyrus__H3K36me3",
    "Brain_Angular_Gyrus__H3K4me1", "Brain_Angular_Gyrus__H3K4me3",
    "Brain_Angular_Gyrus__H3K9ac",
    "Brain_Anterior_Caudate__H3K27ac", "Brain_Anterior_Caudate__H3K36me3",
    "Brain_Anterior_Caudate__H3K4me1", "Brain_Anterior_Caudate__H3K4me3",
    "Brain_Anterior_Caudate__H3K9ac",
    "Brain_Cingulate_Gyrus__H3K27ac", "Brain_Cingulate_Gyrus__H3K36me3",
    "Brain_Cingulate_Gyrus__H3K4me1", "Brain_Cingulate_Gyrus__H3K4me3",
    "Brain_Cingulate_Gyrus__H3K9ac",
    "Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac",
    "Brain_Dorsolateral_Prefrontal_Cortex__H3K36me3",
    "Brain_Dorsolateral_Prefrontal_Cortex__H3K4me1",
    "Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3",
    "Brain_Dorsolateral_Prefrontal_Cortex__H3K9ac",
    "Brain_Germinal_Matrix__H3K36me3",
    "Brain_Germinal_Matrix__H3K4me1",
    "Brain_Germinal_Matrix__H3K4me3",
    "Brain_Hippocampus_Middle__H3K27ac",
    "Brain_Hippocampus_Middle__H3K36me3",
    "Brain_Hippocampus_Middle__H3K4me1",
    "Brain_Hippocampus_Middle__H3K4me3",
    "Brain_Inferior_Temporal_Lobe__H3K27ac",
    "Brain_Inferior_Temporal_Lobe__H3K36me3",
    "Brain_Inferior_Temporal_Lobe__H3K4me1",
    "Brain_Inferior_Temporal_Lobe__H3K4me3",
    "Brain_Inferior_Temporal_Lobe__H3K9ac",
    "Brain_Substantia_Nigra__H3K27ac",
    "Brain_Substantia_Nigra__H3K36me3",
    "Brain_Substantia_Nigra__H3K4me1",
    "Brain_Substantia_Nigra__H3K4me3",
    "Brain_Substantia_Nigra__H3K9ac",
    "Fetal_Adrenal_Gland__DNase", "Fetal_Adrenal_Gland__H3K27ac",
    "Fetal_Adrenal_Gland__H3K36me3", "Fetal_Adrenal_Gland__H3K4me1",
    "Fetal_Adrenal_Gland__H3K4me3",
    "Fetal_Brain_Female__DNase", "Fetal_Brain_Female__H3K36me3",
    "Fetal_Brain_Female__H3K4me1", "Fetal_Brain_Female__H3K4me3",
    "Fetal_Brain_Male__DNase", "Fetal_Brain_Male__H3K36me3",
    "Fetal_Brain_Male__H3K4me1", "Fetal_Brain_Male__H3K4me3",
    "Foreskin_Fibroblast_Primary_Cells_skin01__H3K27ac",
    "HMEC_Mammary_Epithelial_Primary_Cells__H3K27ac",
    "Ovary__DNase", "Ovary__H3K27ac", "Ovary__H3K36me3",
    "Ovary__H3K4me1", "Ovary__H3K4me3",
    "Primary_mononuclear_cells_from_peripheral_blood__H3K4me1",
    "Rectal_Mucosa_Donor_29__H3K4me1"
  ),
  source_prefix = c(
    "Roadmap.224.",
    "Roadmap.60.", "Roadmap.139.", "Roadmap.227.", "Roadmap.315.",
    "Roadmap.373.",
    "Roadmap.61.", "Roadmap.140.", "Roadmap.228.", "Roadmap.316.",
    "Roadmap.374.",
    "Roadmap.62.", "Roadmap.141.", "Roadmap.229.", "Roadmap.317.",
    "Roadmap.375.",
    "Roadmap.65.", "Roadmap.145.", "Roadmap.233.", "Roadmap.321.",
    "Roadmap.377.",
    "Roadmap.142.", "Roadmap.230.", "Roadmap.318.",
    "Roadmap.63.", "Roadmap.143.", "Roadmap.231.", "Roadmap.319.",
    "Roadmap.64.", "Roadmap.144.", "Roadmap.232.", "Roadmap.320.",
    "Roadmap.376.",
    "Roadmap.66.", "Roadmap.146.", "Roadmap.234.", "Roadmap.322.",
    "Roadmap.378.",
    "Roadmap.13.", "Roadmap.71.", "Roadmap.152.", "Roadmap.240.",
    "Roadmap.328.",
    "Roadmap.15.", "Roadmap.154.", "Roadmap.242.", "Roadmap.330.",
    "Roadmap.14.", "Roadmap.153.", "Roadmap.241.", "Roadmap.329.",
    "Roadmap.51.", "Roadmap.98.",
    "Roadmap.27.", "Roadmap.83.", "Roadmap.169.", "Roadmap.257.",
    "Roadmap.345.",
    "Roadmap.223.",
    "Roadmap.261."
  )
)
chromatin_specs[, prefix := file.path(chrom_dir, source_prefix)]

custom_headers <- c(
  "PI", "ASC1", "ASC2", "END", "exCA1", "exCA3", "exDG", "exPFC1",
  "exPFC2", "GABA1", "GABA2", "MG", "NSC", "ODC", "OPC",
  "PI_x_ASC1", "PI_x_ASC2", "PI_x_END", "PI_x_exCA1", "PI_x_exCA3",
  "PI_x_exDG", "PI_x_exPFC1", "PI_x_exPFC2", "PI_x_GABA1", "PI_x_GABA2",
  "PI_x_MG", "PI_x_NSC", "PI_x_ODC", "PI_x_OPC"
)
custom_published_labels <- c(
  "PI Genes", "ASC1", "ASC2", "Endothelial Cells", "exCA1", "exCA3",
  "exDG", "exPFC1", "exPFC2", "GABA1", "GABA2", "Microglia",
  "Neuronal Stem Cells", "Oligodendrocytes", "Oligodendrocyte Precursor Cells",
  "PI x ASC1", "PI x ASC2", "PI x Endothelial Cells", "PI x exCA1",
  "PI x exCA3", "PI x exDG", "PI x exPFC1", "PI x exPFC2", "PI x GABA1",
  "PI x GABA2", "PI x Microglia", "PI x Neuronal Stem Cells",
  "PI x Oligodendrocytes", "PI x Oligodendrocyte Precursor Cells"
)
custom_specs <- data.table(
  annotation_class = "custom",
  published_label = custom_published_labels,
  source_label = custom_headers,
  source_prefix = custom_headers,
  prefix = custom_prefix
)

external_specs <- rbind(expression_specs, chromatin_specs)
all_specs <- rbind(
  external_specs[, .(annotation_class, published_label, source_label, source_prefix, prefix)],
  custom_specs
)

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
  "F_COMP~~SCZ_PGC", "F_INT~~SCZ_PGC", "F_SUD~~SCZ_PGC",
  "F_COMP~~SCZ_FG",  "F_INT~~SCZ_FG",  "F_SUD~~SCZ_FG"
)

prefix_qc <- function(prefix, require_annot = TRUE) {
  ld <- Sys.glob(paste0(prefix, "*l2.ldscore*"))
  m <- Sys.glob(paste0(prefix, "*l2.M_5_50"))
  an <- Sys.glob(paste0(prefix, "*annot.gz"))
  data.table(
    prefix = prefix,
    n_ldscore = length(ld), n_M_5_50 = length(m), n_annot = length(an),
    ldscore_22 = length(ld) == 22L, M_5_50_22 = length(m) == 22L,
    annot_22 = if (require_annot) length(an) == 22L else NA
  )
}

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
if (!all(header_qc$exists & header_qc$has_SNP & header_qc$has_A1 & header_qc$has_A2 &
         header_qc$has_Z & header_qc$has_N)) {
  stop("One or more D6B summary-statistic inputs failed existence/header QC")
}

reference_qc <- rbind(
  data.table(reference = "BaselineLD_v2.2", prefix = baseline_dir,
             n_ldscore = sum(grepl("baselineLD\\.[0-9]+\\.l2\\.ldscore\\.gz$", list.files(baseline_dir))),
             n_M_5_50 = sum(grepl("baselineLD\\.[0-9]+\\.l2\\.M_5_50$", list.files(baseline_dir))),
             n_annot = sum(grepl("baselineLD\\.[0-9]+\\.annot\\.gz$", list.files(baseline_dir)))),
  data.table(reference = "weights_hm3_no_MHC", prefix = wld,
             n_ldscore = sum(grepl("weights\\.hm3_noMHC\\.[0-9]+\\.l2\\.ldscore\\.gz$", list.files(wld))),
             n_M_5_50 = NA_integer_, n_annot = NA_integer_),
  data.table(reference = "frq", prefix = frq,
             n_ldscore = NA_integer_, n_M_5_50 = NA_integer_,
             n_annot = sum(grepl("1000G\\.EUR\\.QC\\.[0-9]+\\.frq$", list.files(frq))))
)
fwrite(reference_qc, file.path(out_dir, "reference_qc.tsv"), sep = "\t")

spec_qc <- rbindlist(lapply(all_specs$prefix, prefix_qc))
all_specs[, c("n_ldscore", "n_M_5_50", "n_annot", "ldscore_22", "M_5_50_22", "annot_22") :=
  spec_qc[, .(n_ldscore, n_M_5_50, n_annot, ldscore_22, M_5_50_22, annot_22)]]
fwrite(all_specs, file.path(out_dir, "D6B_annotation_manifest.tsv"), sep = "\t")

if (!all(reference_qc[reference == "BaselineLD_v2.2", n_ldscore] == 22L,
        reference_qc[reference == "BaselineLD_v2.2", n_M_5_50] == 22L,
        reference_qc[reference == "BaselineLD_v2.2", n_annot] == 22L,
        reference_qc[reference == "weights_hm3_no_MHC", n_ldscore] == 22L,
        reference_qc[reference == "frq", n_annot] == 22L,
        all(all_specs$ldscore_22), all(all_specs$M_5_50_22), all(all_specs$annot_22))) {
  stop("D6B baseline, external, custom, weight, or frequency references are incomplete")
}

if (identical(Sys.getenv("D6B_PREFLIGHT"), "1")) {
  cat(paste0("D6B_PREFLIGHT_OK n_external=", nrow(external_specs),
             " n_custom=", nrow(custom_specs),
             " n_total_extra_prefixes=", length(c(external_specs$prefix, custom_prefix)),
             "\n"), flush = TRUE)
  quit(save = "no", status = 0)
}

config <- c(
  paste0("generated_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("GenomicSEM_version=", as.character(packageVersion("GenomicSEM"))),
  paste0("R_version=", R.version.string),
  "analysis=D6B published PI x brain-cell Stratified Genomic SEM",
  "protocol_source=Grotzinger_2022_Nature_Genetics_Supplementary_Table_10",
  paste0("traits=", paste(names(traits), collapse = ",")),
  "baseline=97 BaselineLD v2.2 annotations",
  "selected_expression=24 (GTEx/Franke; brain/endocrine plus published controls)",
  "selected_chromatin=60 (Roadmap; published tissue/mark selection plus published controls)",
  "custom=29 (PI, 14 cell classes, 14 PI-by-cell intersections)",
  "conditioning_annotation_total=210 (97 + 24 + 60 + 29)",
  "published_binary_enrichment_family=168 (55 BaselineLD binary + 24 expression + 60 chromatin + 29 custom)",
  "baseline_continuous_and_flanking=42 retained for genome-wide conditioning and excluded from enrichment output by official settings",
  "rectal_mucosa_mapping=published label is donor-unspecified; deterministic source mapping uses Roadmap.261 (Rectal Mucosa Donor 29); Roadmap.262 is retained in the audit as the unresolved alternative",
  paste0("external_prefix_count=", nrow(external_specs)),
  paste0("custom_prefix=", custom_prefix),
  paste0("parameters=", paste(params, collapse = ",")),
  "sample_prev=NA for all traits; effective-N sumstats are used",
  "n.blocks=200",
  "exclude_cont=TRUE",
  "tau=FALSE",
  "fix=regressions",
  "std.lv=FALSE",
  "base=TRUE",
  "rm_flank=TRUE",
  "boundary=resid_OCD > 0.0001; inherited from D5i; no new boundary",
  "inference_source=LDSC 200-block jackknife V and GenomicSEM enrich sandwich covariance",
  "additional_sem_block_refit_gate=none",
  "published_alpha_168=0.05/168",
  "covariance_protocol_alpha_155=0.05/155",
  "multiple_testing=report both published 168-family and 155-surviving covariance-protocol thresholds; no custom annotation-specific gate"
)
writeLines(config, file.path(out_dir, "D6B_CONFIG.txt"))
writeLines(model, file.path(out_dir, "D6B_model_syntax.txt"))
fwrite(data.table(trait = names(traits), path = unname(traits)), file.path(out_dir, "trait_manifest.tsv"), sep = "\t")
fwrite(data.table(parameter = params), file.path(out_dir, "target_parameter_manifest.tsv"), sep = "\t")

fail <- function(e) {
  writeLines(c(config, "status=ERROR", paste0("error=", conditionMessage(e)),
               paste0("finished_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))), status_file)
  stop(e)
}

resume_enrich <- identical(Sys.getenv("D6B_RESUME_ENRICH"), "1") &&
  file.exists(file.path(out_dir, "D6B_S_LDSCOutput.rds")) &&
  file.exists(file.path(out_dir, "D6B_enrich_raw.rds"))

tryCatch({
  if (resume_enrich) {
    cat("D6B_RESUME_ENRICH_START\n", flush = TRUE)
    S_LDSCOutput <- readRDS(file.path(out_dir, "D6B_S_LDSCOutput.rds"))
    E <- readRDS(file.path(out_dir, "D6B_enrich_raw.rds"))
    if (is.null(S_LDSCOutput$S) || is.null(S_LDSCOutput$V)) stop("saved s_ldsc object has no S/V")
    n_s <- length(S_LDSCOutput$S)
    cat(paste0("D6B_RESUME_ENRICH_LOADED n_matrices=", n_s, "\n"), flush = TRUE)
  } else {
    cat("D6B_S_LDSC_START\n", flush = TRUE)
    extra_ld <- c(external_specs$prefix, custom_prefix)
    S_LDSCOutput <- s_ldsc(
      traits = unname(traits),
      sample.prev = rep(NA_real_, length(traits)),
      population.prev = rep(NA_real_, length(traits)),
      ld = c(baseline_dir, extra_ld),
      wld = wld,
      frq = frq,
      trait.names = names(traits),
      n.blocks = 200L,
      ldsc.log = file.path(out_dir, "D6B_s_ldsc"),
      exclude_cont = TRUE
    )
    if (is.null(S_LDSCOutput$S) || is.null(S_LDSCOutput$V)) stop("s_ldsc returned no S/V")
    saveRDS(S_LDSCOutput, file.path(out_dir, "D6B_S_LDSCOutput.rds"), compress = FALSE)
    n_s <- length(S_LDSCOutput$S)
    cat(paste0("D6B_S_LDSC_DONE n_matrices=", n_s, "\n"), flush = TRUE)

    cat("D6B_ENRICH_START\n", flush = TRUE)
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
    saveRDS(E, file.path(out_dir, "D6B_enrich_raw.rds"), compress = FALSE)
    cat("D6B_ENRICH_DONE\n", flush = TRUE)
  }

  pieces <- lapply(seq_along(E), function(i) {
    z <- as.data.table(E[[i]])
    z[, target_parameter := if (i <= length(params)) params[[i]] else "BASELINE"]
    z
  })
  enrichment <- rbindlist(pieces, fill = TRUE)
  if (!"Annotation" %in% names(enrichment)) stop("enrich returned no Annotation column")
  enrichment[, annotation_raw := as.character(Annotation)]
  enrichment[, annotation_key := sub("L2$", "", annotation_raw)]
  enrichment[, annotation_basename := basename(sub("/$", "", annotation_raw))]
  enrichment[, annotation_class := "unclassified"]
  enrichment[annotation_key %chin% custom_headers, annotation_class := "custom"]
  for (i in seq_len(nrow(external_specs))) {
    bn <- basename(sub("/$", "", external_specs$prefix[i]))
    enrichment[annotation_basename == bn | annotation_raw == external_specs$prefix[i],
               annotation_class := external_specs$annotation_class[i]]
  }
  enrichment[annotation_class == "unclassified" & grepl("L2$", annotation_raw),
             annotation_class := "baseline"]
  enrichment[target_parameter == "BASELINE", annotation_class := "baseline"]
  enrichment[, published_label := NA_character_]
  for (i in seq_len(nrow(all_specs))) {
    bn <- basename(sub("/$", "", all_specs$prefix[i]))
    enrichment[annotation_key == all_specs$source_label[i] |
                 annotation_raw == all_specs$source_prefix[i] |
                 annotation_basename == bn |
                 annotation_raw == all_specs$prefix[i],
               published_label := all_specs$published_label[i]]
  }
  enrichment[target_parameter == "BASELINE" & is.na(published_label), published_label := "BaselineLD"]
  enrichment[, is_nonbase := target_parameter != "BASELINE" & !is.na(Annotation)]
  enrichment[, alpha_168 := 0.05 / 168]
  enrichment[, alpha_155 := 0.05 / 155]
  enrichment[, Bonferroni_significant_168 := is_nonbase & !is.na(Enrichment_p_value) & Enrichment_p_value < alpha_168]
  enrichment[, Bonferroni_significant_155 := is_nonbase & !is.na(Enrichment_p_value) & Enrichment_p_value < alpha_155]
  enrichment[, Annotation_model_warning := !is.na(Warning) & trimws(as.character(Warning)) != "0"]
  enrichment[, Annotation_model_error := !is.na(Error) & trimws(as.character(Error)) != "0"]
  fwrite(enrichment, file.path(out_dir, "D6B_enrichment_long.tsv"), sep = "\t", na = "NA")

  tested <- enrichment[is_nonbase == TRUE]
  n_annotation_tests_actual <- length(unique(tested$Annotation))
  sig168 <- tested[Bonferroni_significant_168 == TRUE]
  sig155 <- tested[Bonferroni_significant_155 == TRUE]
  fwrite(sig168, file.path(out_dir, "D6B_enrichment_bonferroni_168_significant.tsv"), sep = "\t", na = "NA")
  fwrite(sig155, file.path(out_dir, "D6B_enrichment_bonferroni_155_significant.tsv"), sep = "\t", na = "NA")
  fwrite(tested[annotation_class == "custom"], file.path(out_dir, "D6B_custom_enrichment.tsv"), sep = "\t", na = "NA")

  parameter_summary <- tested[, .(
    n_rows = .N,
    n_finite_enrichment = sum(is.finite(Enrichment)),
    n_finite_p = sum(is.finite(Enrichment_p_value)),
    n_bonferroni_significant_168 = sum(Bonferroni_significant_168, na.rm = TRUE),
    n_bonferroni_significant_155 = sum(Bonferroni_significant_155, na.rm = TRUE),
    n_annotation_warnings = sum(Annotation_model_warning, na.rm = TRUE),
    n_annotation_errors = sum(Annotation_model_error, na.rm = TRUE),
    min_p = suppressWarnings(min(Enrichment_p_value, na.rm = TRUE))
  ), by = target_parameter]
  fwrite(parameter_summary, file.path(out_dir, "D6B_parameter_summary.tsv"), sep = "\t", na = "NA")

  class_summary <- tested[, .(
    n_rows = .N,
    n_unique_annotations = uniqueN(Annotation),
    n_finite_p = sum(is.finite(Enrichment_p_value)),
    n_sig_168 = sum(Bonferroni_significant_168, na.rm = TRUE),
    n_sig_155 = sum(Bonferroni_significant_155, na.rm = TRUE)
  ), by = annotation_class]
  fwrite(class_summary, file.path(out_dir, "D6B_annotation_class_summary.tsv"), sep = "\t", na = "NA")

  fwrite(data.table(
    n_s_ldsc_matrices = n_s,
    n_enrich_list_elements = length(E),
    n_rows_nonbase = nrow(tested),
    n_annotation_tests_actual = n_annotation_tests_actual,
    n_annotation_tests_published = 168L,
    n_parameters = length(params),
    alpha_168 = 0.05 / 168,
    alpha_155 = 0.05 / 155,
    n_significant_168 = nrow(sig168),
    n_significant_155 = nrow(sig155),
    n_warnings = sum(tested$Annotation_model_warning, na.rm = TRUE),
    n_errors = sum(tested$Annotation_model_error, na.rm = TRUE)
  ), file.path(out_dir, "D6B_family_summary.tsv"), sep = "\t", na = "NA")

  writeLines(c(config,
               "status=COMPLETED",
               paste0("n_s_ldsc_matrices=", n_s),
               paste0("n_enrich_list_elements=", length(E)),
               paste0("n_rows_nonbase=", nrow(tested)),
               paste0("n_annotation_tests_actual=", n_annotation_tests_actual),
               paste0("n_significant_168=", nrow(sig168)),
               paste0("n_significant_155=", nrow(sig155)),
               paste0("finished_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))), status_file)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  cat("D6B_COMPLETE\n", flush = TRUE)
}, error = fail)

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) args[[1]] else "."

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

result_root <- if (dir.exists(file.path(base, "results", "ascertainment_pilot", "p25"))) {
  file.path(base, "results", "ascertainment_pilot")
} else {
  file.path(base, "results")
}

out_dir <- file.path(result_root, "p9_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

elig <- read_tsv(file.path(result_root, "p25", "gates",
                           "family_common_eligibility.tsv"))
elig$all_rg2_lt_0.5 <- tolower(as.character(elig$all_rg2_lt_0.5)) == "true"
prov <- read_tsv(file.path(base, "manuscript", "supplementary_tables",
                           "Supplementary_Table_2_psychiatric_auxiliary_provenance_and_h2.tsv"))

elig$rg2_threshold <- 0.5
elig$max_rg2 <- elig$max_abs_rg^2
elig$h2_z <- prov$z_h2[match(elig$auxiliary, prov$trait)]
elig$h2_eligible_z_gt_6 <- elig$h2_z > 6
elig$selected <- elig$h2_eligible_z_gt_6 & elig$all_rg2_lt_0.5
elig$exclusion_reason <- ifelse(
  !elig$h2_eligible_z_gt_6,
  "heritability_z_not_greater_than_6",
  ifelse(!elig$all_rg2_lt_0.5,
         "at_least_one_target_auxiliary_rg_squared_not_below_0.5",
         "included")
)
elig$source_identifier <- prov$source_identifier[match(elig$auxiliary, prov$trait)]
elig$source_status <- prov$source_status[match(elig$auxiliary, prov$trait)]

psy <- data.frame(
  analysis_family = elig$family,
  candidate_key = elig$auxiliary,
  candidate_trait = elig$auxiliary,
  definitions_tested = elig$definitions_tested,
  h2_z = elig$h2_z,
  h2_eligible_z_gt_6 = elig$h2_eligible_z_gt_6,
  max_abs_target_aux_rg = elig$max_abs_rg,
  max_target_aux_rg2 = elig$max_rg2,
  rg2_threshold = elig$rg2_threshold,
  all_target_aux_rg2_below_threshold = elig$all_rg2_lt_0.5,
  selected = elig$selected,
  exclusion_reason = elig$exclusion_reason,
  source_identifier = elig$source_identifier,
  source_status = elig$source_status,
  stringsAsFactors = FALSE
)

neu_manifest <- read_tsv(file.path(result_root, "p5",
                                   "neurological_auxiliary_candidate_manifest.tsv"))
neu_gate <- read_tsv(file.path(result_root, "p5",
                               "neurological_auxiliary_h2_gate.tsv"))
neu <- merge(neu_manifest, neu_gate, by = c("key", "trait"), sort = FALSE)
neu <- neu[match(neu_manifest$key, neu$key), ]
neu_selected <- neu$gate == "PASS"

neu_out <- data.frame(
  analysis_family = "NEUROLOGICAL",
  candidate_key = neu$key,
  candidate_trait = neu$trait,
  definitions_tested = NA_integer_,
  h2_z = neu$z_h2,
  h2_eligible_z_gt_6 = neu$z_h2 > 6,
  max_abs_target_aux_rg = NA_real_,
  max_target_aux_rg2 = NA_real_,
  rg2_threshold = NA_real_,
  all_target_aux_rg2_below_threshold = NA,
  selected = neu_selected,
  exclusion_reason = ifelse(neu_selected, "included", "heritability_z_not_greater_than_6"),
  source_identifier = paste(neu$accession, neu$source_url, sep = "; "),
  source_status = ifelse(neu$public_summary_statistics == "yes", "IDENTIFIED", "UNAVAILABLE"),
  stringsAsFactors = FALSE
)

audit <- rbind(psy, neu_out)
if (anyNA(audit$h2_z)) stop("Missing h2 Z value in candidate audit")
if (any(audit$selected & audit$exclusion_reason != "included")) {
  stop("Selected candidate has an exclusion reason")
}

write.table(
  audit,
  file.path(out_dir, "auxiliary_candidate_selection_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

summary <- aggregate(
  selected ~ analysis_family,
  audit,
  function(x) c(candidate_n = length(x), selected_n = sum(x), excluded_n = sum(!x))
)
summary <- data.frame(
  analysis_family = summary$analysis_family,
  candidate_n = summary$selected[, "candidate_n"],
  selected_n = summary$selected[, "selected_n"],
  excluded_n = summary$selected[, "excluded_n"],
  stringsAsFactors = FALSE
)
write.table(
  summary,
  file.path(out_dir, "auxiliary_candidate_selection_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Wrote", nrow(audit), "candidate-family rows\n")

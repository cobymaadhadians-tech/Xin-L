#!/usr/bin/env Rscript

root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
out_dir <- file.path(root, "results/p9_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sources <- list(
  psychiatric_target_target = file.path(root, "results/p2/ldsc_rg"),
  psychiatric_target_auxiliary = file.path(root, "results/p25/target_aux_rg"),
  psychiatric_auxiliary_auxiliary = file.path(root, "results/p3/aux_rg"),
  neurological = file.path(root, "results/p5/neurological_auxiliary/ldsc_rg")
)

extract_number <- function(pattern, text, field, path) {
  hit <- regexec(pattern, text, perl = TRUE)
  value <- regmatches(text, hit)[[1L]]
  if (length(value) < 2L) stop("could not parse ", field, " from ", path)
  as.numeric(value[[2L]])
}

parse_log <- function(path, scope) {
  lines <- readLines(path, warn = FALSE)
  text <- paste(lines, collapse = "\n")
  name <- sub("\\.log$", "", basename(path))
  traits <- strsplit(name, "__", fixed = TRUE)[[1L]]
  if (length(traits) != 2L) stop("unexpected pair filename: ", path)
  section_start <- regexpr("Genetic Covariance\\n-+", text, perl = TRUE)
  if (section_start < 0L) stop("missing genetic covariance section: ", path)
  covariance_text <- substring(text, section_start)
  next_section <- regexpr("\\nGenetic Correlation\\n-+", covariance_text, perl = TRUE)
  if (next_section < 0L) stop("missing genetic correlation section: ", path)
  covariance_text <- substring(covariance_text, 1L, next_section - 1L)
  intercept_match <- regexec(
    "Intercept:[[:space:]]*([-+0-9.eE]+)[[:space:]]*\\(([-+0-9.eE]+)\\)",
    covariance_text, perl = TRUE
  )
  intercept_values <- regmatches(covariance_text, intercept_match)[[1L]]
  if (length(intercept_values) != 3L) stop("could not parse genetic covariance intercept: ", path)
  constrained_flag <- grepl("--intercept-gencov|--no-intercept", text, perl = TRUE)
  data.frame(
    scope = scope,
    trait_1 = traits[[1L]],
    trait_2 = traits[[2L]],
    genetic_covariance_intercept = as.numeric(intercept_values[[2L]]),
    genetic_covariance_intercept_se = as.numeric(intercept_values[[3L]]),
    genetic_correlation = extract_number(
      "Genetic Correlation:[[:space:]]*([-+0-9.eE]+)", text, "genetic correlation", path
    ),
    genetic_correlation_se = extract_number(
      "Genetic Correlation:[[:space:]]*[-+0-9.eE]+[[:space:]]*\\(([-+0-9.eE]+)\\)",
      text, "genetic correlation SE", path
    ),
    cross_trait_intercept_mode = if (constrained_flag) "CONSTRAINED_FLAG_PRESENT" else "FREE_ESTIMATED",
    ldsc_version = if (grepl("Version 1.0.1", text, fixed = TRUE)) "1.0.1" else NA_character_,
    log_path = sub(paste0("^", root, "/"), "", path),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (scope in names(sources)) {
  paths <- sort(list.files(sources[[scope]], pattern = "\\.log$", full.names = TRUE))
  if (!length(paths)) stop("no LDSC logs found for ", scope)
  for (path in paths) rows[[length(rows) + 1L]] <- parse_log(path, scope)
}
result <- do.call(rbind, rows)
if (anyDuplicated(result[c("scope", "trait_1", "trait_2")])) stop("duplicate LDSC pair rows")
write.table(
  result, file.path(out_dir, "ldsc_cross_trait_intercepts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("P9_LDSC_INTERCEPT_EXTRACTION_COMPLETE\n")

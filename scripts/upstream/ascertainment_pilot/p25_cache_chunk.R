#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("target key, chunk id, and chunk size required")
target_key <- args[[1L]]
chunk_id <- as.integer(args[[2L]])
chunk_size <- as.integer(args[[3L]])
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))
if (!target_key %in% names(p25_panels)) stop("unsupported target: ", target_key)

n_block <- 200L
aux <- p25_panels[[target_key]]
traits <- c(target_key, aux)
start_block <- chunk_id * chunk_size + 1L
end_block <- min(n_block, (chunk_id + 1L) * chunk_size)
if (start_block > n_block) stop("chunk starts after final block")
block_ids <- start_block:end_block

all_paths <- p25_sumstats_paths(root)
paths <- all_paths[traits]
if (anyNA(paths) || !all(file.exists(paths))) stop("missing input for ", target_key)
munged <- lapply(paths, function(path) fread(path, showProgress = FALSE))
names(munged) <- traits

hmp3 <- Sys.getenv("HM3_SNPLIST", unset = file.path(root, "reference/ldsc/w_hm3.snplist"))
ld_path <- Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr"))
wld_path <- Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC"))
hmp3_snp <- fread(hmp3, showProgress = FALSE)$SNP
blocks <- split(hmp3_snp, cut(seq_along(hmp3_snp), breaks = n_block, labels = FALSE))

rg_chunk <- array(NA_real_, dim = c(length(traits), length(traits), length(block_ids)), dimnames = list(traits, traits, paste0("block_", block_ids)))
h2_chunk <- matrix(NA_real_, nrow = length(block_ids), ncol = length(traits), dimnames = list(paste0("block_", block_ids), traits))
gcov_chunk <- array(NA_real_, dim = c(length(traits), length(traits), length(block_ids)), dimnames = list(traits, traits, paste0("block_", block_ids)))

for (j in seq_along(block_ids)) {
  block <- block_ids[[j]]
  remain <- lapply(munged, function(gwas) gwas[!is.na(N) & !SNP %in% blocks[[block]]])
  names(remain) <- paste0("GWAS_", traits, "_remainblock")
  rg_res <- pleioh2g:::ldsc_rg(
    munged_sumstats = remain,
    sample_prev = NA,
    population_prev = NA,
    ld = ld_path,
    wld = wld_path,
    n_blocks = 200,
    chisq_max = NA,
    chr_filter = 1:22
  )
  rg <- matrix(0, nrow = length(traits), ncol = length(traits), dimnames = list(traits, traits))
  lower <- which(lower.tri(rg, diag = FALSE), arr.ind = TRUE)
  rg[lower] <- rg_res$rg$rg
  rg <- rg + t(rg)
  diag(rg) <- 1
  rg_chunk[, , j] <- rg
  h2_chunk[j, ] <- rg_res$h2$h2_observed
  gcov_chunk[, , j] <- rg_res$raw$S
}

out_dir <- file.path(root, "results/p25/phbc_cache", target_key)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("chunk_%03d.rds", chunk_id))
part_file <- paste0(out_file, ".part")
saveRDS(list(
  target = target_key,
  traits = traits,
  auxiliaries = aux,
  n_block = n_block,
  block_ids = block_ids,
  rgarray = rg_chunk,
  h2array = h2_chunk,
  gcovarray = gcov_chunk,
  pleioh2g_version = as.character(packageVersion("pleioh2g"))
), part_file)
file.rename(part_file, out_file)
cat("P25_CACHE_CHUNK_COMPLETE\t", target_key, "\t", chunk_id, "\t", start_block, "-", end_block, "\n", sep = "")

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("chunk id required")
chunk_id <- as.integer(args[[1]])
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pleioh2g)
})
source(file.path(root, "scripts/upstream/ascertainment_pilot/p25_panel_config.R"))

target <- "SCZ_PGC2022"
aux <- p25_panels[[target]]
traits <- c(target, aux)
n_block <- 200L
chunk_size <- 10L
block_ids <- (chunk_id * chunk_size + 1L):((chunk_id + 1L) * chunk_size)
if (chunk_id < 0L || chunk_id > 19L) stop("chunk must be 0..19")

paths <- p25_sumstats_paths(root)[traits]
paths[[target]] <- file.path(root, "results/p11_scz_neff_audit/input/SCZ_PGC2022_Nx2.sumstats.gz")
if (anyNA(paths) || !all(file.exists(paths))) stop("missing corrected cache input")
munged <- lapply(paths, fread, showProgress = FALSE)
names(munged) <- traits

hmp3_snp <- fread(Sys.getenv("HM3_SNPLIST", unset = file.path(root, "reference/ldsc/w_hm3.snplist")), showProgress = FALSE)$SNP
blocks <- split(hmp3_snp, cut(seq_along(hmp3_snp), breaks = n_block, labels = FALSE))
rg_chunk <- array(NA_real_, c(length(traits), length(traits), length(block_ids)),
                  dimnames = list(traits, traits, paste0("block_", block_ids)))
h2_chunk <- matrix(NA_real_, length(block_ids), length(traits),
                   dimnames = list(paste0("block_", block_ids), traits))
gcov_chunk <- array(NA_real_, c(length(traits), length(traits), length(block_ids)),
                    dimnames = list(traits, traits, paste0("block_", block_ids)))

for (j in seq_along(block_ids)) {
  b <- block_ids[[j]]
  remain <- lapply(munged, function(gwas) gwas[!is.na(N) & !SNP %in% blocks[[b]]])
  names(remain) <- paste0("GWAS_", traits, "_remainblock")
  z <- pleioh2g:::ldsc_rg(
    munged_sumstats = remain, sample_prev = NA, population_prev = NA,
    ld = Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr")),
    wld = Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC")),
    n_blocks = 200, chisq_max = NA, chr_filter = 1:22
  )
  rg <- matrix(0, length(traits), length(traits), dimnames = list(traits, traits))
  lower <- which(lower.tri(rg), arr.ind = TRUE)
  rg[lower] <- z$rg$rg
  rg <- rg + t(rg)
  diag(rg) <- 1
  rg_chunk[, , j] <- rg
  h2_chunk[j, ] <- z$h2$h2_observed
  gcov_chunk[, , j] <- z$raw$S
}

out <- file.path(root, "results/p11_scz_neff_audit/phbc_cache")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(target = target, traits = traits, auxiliaries = aux, n_block = n_block,
             block_ids = block_ids, rgarray = rg_chunk, h2array = h2_chunk,
             gcovarray = gcov_chunk, pleioh2g_version = as.character(packageVersion("pleioh2g")),
             neff_convention = "2x_source_NEFF"),
        file.path(out, sprintf("chunk_%03d.rds", chunk_id)))

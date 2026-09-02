#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("chunk id and chunk size required")
chunk_id <- as.integer(args[[1L]])
chunk_size <- as.integer(args[[2L]])
root <- normalizePath(Sys.getenv("ANALYSIS_ROOT", unset = getwd()), mustWork = FALSE)
suppressPackageStartupMessages({library(data.table); library(pleioh2g)})
target <- "EPILEPSY_FINNGEN_R13"
aux <- c("NEURO_PD_GCST009324", "NEURO_ALS_GCST90027164", "NEURO_MIGRAINE_GCST90271641",
         "NEURO_ISCHEMIC_STROKE_GCST90104540", "NEURO_RLS_GCST90435387")
traits <- c(target, aux)
n_block <- 200L
start_block <- chunk_id * chunk_size + 1L
end_block <- min(n_block, (chunk_id + 1L) * chunk_size)
block_ids <- start_block:end_block
paths <- c(setNames(file.path(root, "results/p5/ldsc_input", paste0(target, ".sumstats.gz")), target),
           setNames(file.path(root, "results/p5/neurological_auxiliary/ldsc_input", paste0(aux, ".sumstats.gz")), aux))
munged <- lapply(paths, fread, showProgress = FALSE)
names(munged) <- traits
hmp3 <- Sys.getenv("HM3_SNPLIST", unset = file.path(root, "reference/ldsc/w_hm3.snplist"))
ld_path <- Sys.getenv("LDSC_LD_PATH", unset = file.path(root, "reference/ldsc/eur_w_ld_chr"))
wld_path <- Sys.getenv("LDSC_WEIGHTS_PATH", unset = file.path(root, "reference/ldsc/weights_hm3_noMHC"))
blocks <- split(fread(hmp3, showProgress = FALSE)$SNP,
                cut(seq_len(nrow(fread(hmp3, showProgress = FALSE))), breaks = n_block, labels = FALSE))
rg_chunk <- array(NA_real_, c(length(traits), length(traits), length(block_ids)),
                  dimnames = list(traits, traits, paste0("block_", block_ids)))
h2_chunk <- matrix(NA_real_, length(block_ids), length(traits),
                   dimnames = list(paste0("block_", block_ids), traits))
gcov_chunk <- array(NA_real_, c(length(traits), length(traits), length(block_ids)),
                    dimnames = list(traits, traits, paste0("block_", block_ids)))
for (j in seq_along(block_ids)) {
  block <- block_ids[[j]]
  remain <- lapply(munged, function(gwas) gwas[!is.na(N) & !SNP %in% blocks[[block]]])
  names(remain) <- paste0("GWAS_", traits, "_remainblock")
  x <- pleioh2g:::ldsc_rg(munged_sumstats = remain, sample_prev = NA, population_prev = NA,
                          ld = ld_path, wld = wld_path, n_blocks = 200,
                          chisq_max = NA, chr_filter = 1:22)
  rg <- matrix(0, length(traits), length(traits), dimnames = list(traits, traits))
  lower <- which(lower.tri(rg, diag = FALSE), arr.ind = TRUE)
  rg[lower] <- x$rg$rg
  rg <- rg + t(rg); diag(rg) <- 1
  rg_chunk[, , j] <- rg
  h2_chunk[j, ] <- x$h2$h2_observed
  gcov_chunk[, , j] <- x$raw$S
}
out <- file.path(root, "results/p6/epilepsy_finngen/cache")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
file <- file.path(out, sprintf("chunk_%03d.rds", chunk_id))
part <- paste0(file, ".part")
saveRDS(list(target = target, traits = traits, auxiliaries = aux, n_block = n_block,
             block_ids = block_ids, rgarray = rg_chunk, h2array = h2_chunk,
             gcovarray = gcov_chunk), part)
file.rename(part, file)
cat("P6_EPILEPSY_FINNGEN_CACHE_COMPLETE\t", chunk_id, "\n", sep = "")

p26_panels <- list(
  SCZ_FINNGEN_R13_EXMORE = c(
    "MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD", "ASD",
    "OCD_2025", "AN", "AUD", "CUD_2023_EUR"
  )
)

p26_sumstats_paths <- function(root) {
  c(
    SCZ_FINNGEN_R13_EXMORE = file.path(root, "results/p26_broad_scz/ldsc_input/SCZ_FINNGEN_R13_EXMORE.sumstats.gz"),
    MDD_CLIN_PGC2025 = file.path(root, "results/p1/ldsc_input/MDD_CLIN_PGC2025.sumstats.gz"),
    BD_PGC2021 = file.path(root, "results/p1/ldsc_input/BD_PGC2021.sumstats.gz"),
    PTSD = file.path(root, "data/input/ldsc/PTSD.sumstats.gz"),
    ANX = file.path(root, "data/input/ldsc/ANX.sumstats.gz"),
    ADHD = file.path(root, "data/input/ldsc/ADHD.sumstats.gz"),
    ASD = file.path(root, "results/p3/aux_munged/ASD.sumstats.gz"),
    OCD_2025 = file.path(root, "results/p25/aux_munged/OCD_2025.sumstats.gz"),
    AN = file.path(root, "results/p3/aux_munged/AN.sumstats.gz"),
    AUD = file.path(root, "data/input/ldsc/AUD.sumstats.gz"),
    CUD_2023_EUR = file.path(root, "results/p25/aux_munged/CUD_2023_EUR.sumstats.gz")
  )
}

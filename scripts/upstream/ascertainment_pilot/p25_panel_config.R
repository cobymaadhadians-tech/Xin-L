p25_panels <- list(
  MDD_CLIN_PGC2025 = c("BD_PGC2021", "SCZ_PGC2022", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  MDD_EHR_PGC2025 = c("BD_PGC2021", "SCZ_PGC2022", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  MDD_QUEST_PGC2025 = c("BD_PGC2021", "SCZ_PGC2022", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  MDD_FINNGEN_R13 = c("BD_PGC2021", "SCZ_PGC2022", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  BD_CLIN_PGC4 = c("SCZ_PGC2022", "PTSD", "ANX", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  BD_FINNGEN_R13 = c("SCZ_PGC2022", "PTSD", "ANX", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  SCZ_PGC2022 = c("MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR"),
  SCZ_FINNGEN_R13 = c("MDD_CLIN_PGC2025", "BD_PGC2021", "PTSD", "ANX", "ADHD", "ASD", "OCD_2025", "AN", "AUD", "CUD_2023_EUR")
)

p25_sumstats_paths <- function(root) {
  c(
    MDD_CLIN_PGC2025 = file.path(root, "results/p1/ldsc_input/MDD_CLIN_PGC2025.sumstats.gz"),
    MDD_EHR_PGC2025 = file.path(root, "results/p1/ldsc_input/MDD_EHR_PGC2025.sumstats.gz"),
    MDD_QUEST_PGC2025 = file.path(root, "results/p1/ldsc_input/MDD_QUEST_PGC2025.sumstats.gz"),
    MDD_FINNGEN_R13 = file.path(root, "results/p1/ldsc_input/MDD_FINNGEN_R13.sumstats.gz"),
    BD_CLIN_PGC4 = file.path(root, "results/p1/ldsc_input/BD_CLIN_PGC4.sumstats.gz"),
    BD_FINNGEN_R13 = file.path(root, "results/p1/ldsc_input/BD_FINNGEN_R13.sumstats.gz"),
    BD_PGC2021 = file.path(root, "results/p1/ldsc_input/BD_PGC2021.sumstats.gz"),
    SCZ_PGC2022 = file.path(root, "results/p1/ldsc_input/SCZ_PGC2022.sumstats.gz"),
    SCZ_FINNGEN_R13 = file.path(root, "results/p1/ldsc_input/SCZ_FINNGEN_R13.sumstats.gz"),
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

# Input data

The repository excludes GWAS summary-statistic files. Download the public files listed in `data/supplementary_tables/Supplementary_Table_1_target_GWAS_provenance.tsv` and `Supplementary_Table_2_psychiatric_auxiliary_provenance_and_h2.tsv`, then place project-normalized LDSC inputs under `data/input/ldsc/` when running the upstream scripts.

Reference LD resources are also excluded. Set `HM3_SNPLIST`, `LDSC_LD_PATH`, and `LDSC_WEIGHTS_PATH` when running the R scripts, or use the documented defaults under `reference/ldsc/`.

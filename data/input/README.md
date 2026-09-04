# Input data

The repository excludes GWAS summary-statistic files. The expected inputs and their current archive status are listed in `data/input/input_manifest.tsv`. Download the public files listed in `data/supplementary_tables/Supplementary_Table_1_target_GWAS_provenance.tsv` and `Supplementary_Table_2_psychiatric_auxiliary_provenance_and_h2.tsv`, then place project-normalized LDSC inputs under `data/input/ldsc/` when running the upstream scripts. The source and reference records separate input provenance from analysis-process logging.

Reference LD resources are also excluded and are listed in `data/reference/reference_manifest.tsv`. Set `HM3_SNPLIST`, `LDSC_LD_PATH`, and `LDSC_WEIGHTS_PATH` when running the R scripts, or use the documented defaults under `reference/ldsc/`.

# Upstream analysis code

This directory contains the analysis source recovered from the server project. Files are grouped by their role in the analysis path:

1. `core/` standardizes public GWAS files and prepares LDSC-compatible inputs.
2. `ascertainment_pilot/` computes the psychiatric, neurological, PHBC, sensitivity, SCZ model-swap and observed-regression results.
3. `p26_broad_scz/` computes the broad FinnGen SCZ definition sensitivity analysis.
4. `supergnova/` prepares inputs and provides a portable launcher for the pinned SUPERGNOVA source.

Server-specific Slurm headers, logs, private absolute paths, cached RDS objects, compiled Python files and bundled software environments were excluded. Comments that define an input, output or algorithm were retained. Public input locations are controlled by `ANALYSIS_ROOT`, `HM3_SNPLIST`, `LDSC_LD_PATH` and `LDSC_WEIGHTS_PATH`.

The result tables committed under `data/derived/` and `data/supplementary_tables/` are the review-accessible outputs. Re-running the full upstream analysis requires the public GWAS files, LD-score reference files, LDSC, the R packages listed in `environment/R_packages.tsv`, and the SUPERGNOVA reference panel.

# Core analysis scripts

This repository contains the core scripts for the schizophrenia GWAS representation-sensitivity analysis. It is intentionally limited to the principal statistical workflow and its input checks. Raw GWAS files, derived results, figures, supplementary tables, extended upstream workflows and submission documents are supplied separately in the analysis archive.

## Principal workflow

The main analysis uses 995,685 common SNPs and 200 matched deletion blocks. The retained scripts cover:

- common-SNP summary reconstruction and joint genetic R²/Shapley calculations;
- bias diagnostics and PHBC-style correction;
- conditional null calibration for anchored covariance and sharing-profile tests;
- constrained Gaussian tests of equal population joint sharing;
- LDSC command wrappers and input validation.

The principal scripts are `recompute_common_snp_summaries.py`, `compute_r2_bias_diagnostics.py`, `correct_six_sharing.py`, `calibrate_sharing_null.py`, `calibrate_anchored_null.py` and `test_total_sharing_gaussian.py`.

## Running locally

Run the scripts from an analysis root containing the required `results/` inputs. Install the packages listed under `environment/` and provide the GWAS and LD-reference files through the local paths required by the selected workflow. The repository contains no raw GWAS summary statistics, result tables, credentials, private server files or figure-rendering code.

For the complete reconstruction workflow, upstream scripts and derived data, use the accompanying analysis archive supplied with the manuscript.

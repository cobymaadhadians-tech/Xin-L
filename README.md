# Analysis repository

This directory contains the derived result tables, figure source data, figure-rendering code and upstream analysis source for the Translational Psychiatry submission.

## Scope

The repository does not contain GWAS summary-statistic files. Public source URLs, accessions and phenotype representations are recorded in `data/supplementary_tables/Supplementary_Table_1_target_GWAS_provenance.tsv` and the associated auxiliary-trait table.

The repository excludes raw GWAS summary statistics and private server files. The public source records remain in the supplementary provenance tables. Upstream Python and R source recovered from the analysis server is organized under `scripts/upstream/`; server-specific Slurm headers, caches, compiled files and bundled software environments were removed. Pairwise LDSC audit logs are retained under `results/ldsc_audit/logs/`. `data/input/input_manifest.tsv` and `data/reference/reference_manifest.tsv` record the excluded inputs and reference resources with source identifiers, accessions, URLs and checksums; the LDSC per-file checksums are in `data/reference/ldsc_reference_checksums.tsv`.

## Reproduce the figures

From this directory, install the listed Python packages and run:

```bash
python3 scripts/validate_inputs.py
python3 -m pip install -r environment/requirements.txt
python3 scripts/render_manuscript_figures.py
```

The command writes uncompressed 600-d.p.i. TIFF files to `outputs/main/` and `outputs/supplementary/`. The default output directory is separate from the submission package. A different directory can be supplied with `--output-dir`.

## Analysis path

1. Public GWAS files and phenotype representations are identified in the provenance tables.
2. `scripts/upstream/core/` standardizes source files and prepares LDSC-compatible inputs.
3. `scripts/run_ldsc_h2.sh` and `scripts/run_ldsc_rg.sh` run the LDSC h2 and rg commands against an external LDSC installation and LD reference.
4. `scripts/upstream/ascertainment_pilot/` computes PHBC, block-jackknife, leave-one-out, neurological and SCZ sensitivity results.
5. `scripts/upstream/p26_broad_scz/` computes the broad FinnGen SCZ representation sensitivity analysis.
6. `scripts/upstream/ascertainment_pilot/prepare_supergnova_inputs.py` prepares local-covariance inputs. `scripts/upstream/supergnova/run_supergnova.sh` launches an external SUPERGNOVA installation pinned in `environment/external_software.tsv`.
7. The resulting tables under `data/derived/` and `data/supplementary_tables/` are the referee-accessible outputs, including the six-trait SUPERGNOVA summary in Supplementary Table 21 and the LDSC diagnostic and overlap audit in Supplementary Table 22.
8. `scripts/render_manuscript_figures.py` combines the committed derived tables into Figures 1–3 and Supplementary Figures 1–3. The SCZ local-covariance renderer is kept in `scripts/render_figure3_diagnostics.py` and is called by the main renderer.

## Upstream environment

The upstream scripts require R 4.5.0, Python 3.9 for the pinned SUPERGNOVA environment, LDSC 1.0.1 at Git commit `2fdeeb3b44379408794154993dbd6101b8946b7e`, `pleioh2g` 0.1.3, GenomicSEM 0.0.5, `mvtnorm` 1.3-3 and the packages listed in `environment/`. The local LDSC audit rerun used Python 3.14.3 with a documented compatibility copy of that LDSC source. Reference LD resources are excluded from this repository and must be supplied through the environment variables described in `data/input/README.md`. The local submission snapshot is recorded in `environment/submission_snapshot.tsv`; the workspace has no Git metadata, so no remote tag was created.

The output mapping and input paths are listed in `analysis_manifest.tsv`.

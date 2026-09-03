# Analysis repository

This directory contains the derived result tables, figure source data, figure-rendering code and upstream analysis source for the Translational Psychiatry submission.

## Scope

The repository does not contain GWAS summary-statistic files. Public source URLs, accessions and phenotype definitions are recorded in `data/supplementary_tables/Supplementary_Table_1_target_GWAS_provenance.tsv` and the associated auxiliary-trait table.

The repository excludes raw GWAS summary statistics and private server files. The public source records remain in the supplementary provenance tables. Upstream Python and R source recovered from the analysis server is organized under `scripts/upstream/`; server-specific Slurm headers, logs, caches, compiled files and bundled software environments were removed.

## Reproduce the figures

From this directory, install the listed Python packages and run:

```bash
python3 scripts/validate_inputs.py
python3 -m pip install -r environment/requirements.txt
python3 scripts/render_manuscript_figures.py
```

The command writes uncompressed 600-d.p.i. TIFF files to `outputs/main/` and `outputs/supplementary/`. The default output directory is separate from the submission package. A different directory can be supplied with `--output-dir`.

## Analysis path

1. Public GWAS files and phenotype definitions are identified in the provenance tables.
2. `scripts/upstream/core/` standardizes source files and prepares LDSC-compatible inputs.
3. `scripts/run_ldsc_h2.sh` and `scripts/run_ldsc_rg.sh` run the LDSC h2 and rg commands against an external LDSC installation and LD reference.
4. `scripts/upstream/ascertainment_pilot/` computes PHBC, block-jackknife, leave-one-out, neurological and SCZ sensitivity results.
5. `scripts/upstream/p26_broad_scz/` computes the broad FinnGen SCZ definition sensitivity analysis.
6. `scripts/upstream/ascertainment_pilot/prepare_supergnova_inputs.py` prepares local-covariance inputs. `scripts/upstream/supergnova/run_supergnova.sh` launches an external SUPERGNOVA installation pinned in `environment/external_software.tsv`.
7. The resulting tables under `data/derived/` and `data/supplementary_tables/` are the referee-accessible outputs.
8. `scripts/render_manuscript_figures.py` combines the committed derived tables into Figures 1–3 and Supplementary Figures 1–3. The SCZ local-covariance renderer is kept in `scripts/render_figure3_diagnostics.py` and is called by the main renderer.

## Upstream environment

The upstream scripts require R 4.5.0, Python 3.9 for the pinned SUPERGNOVA environment, LDSC, `pleioh2g` 0.1.3, GenomicSEM 0.0.5 and the packages listed in `environment/`. Reference LD resources are excluded from this repository and must be supplied through the environment variables described in `data/input/README.md`.

The output mapping and input paths are listed in `analysis_manifest.tsv`.

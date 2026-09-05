# Analysis repository

This directory contains the derived result tables, figure source data, figure-rendering code and upstream analysis source for the Molecular Psychiatry submission.

## Scope

The repository does not contain GWAS summary-statistic files. Public source URLs, accessions and phenotype representations are recorded in `data/supplementary_tables/Supplementary_Table_1_target_GWAS_provenance.tsv` and the associated auxiliary-trait table.

The repository excludes raw GWAS summary statistics and private server files. The public source records remain in the supplementary provenance tables. Upstream Python and R source recovered from the analysis server is organized under `scripts/upstream/`; server-specific Slurm headers, caches, compiled files and bundled software environments were removed. Pairwise LDSC audit logs are retained under `results/ldsc_audit/logs/`. `data/input/input_manifest.tsv` and `data/reference/reference_manifest.tsv` record the excluded inputs and reference resources with source identifiers, accessions, URLs and checksums; the LDSC per-file checksums are in `data/reference/ldsc_reference_checksums.tsv`.

## Reproduce the figures

From this directory, install the listed Python packages and run:

```bash
python3 scripts/validate_inputs.py
python3 -m pip install -r environment/requirements.txt
python3 scripts/render_manuscript_figures.py
python3 scripts/render_figure4_downstream.py --outdir outputs/supplementary --pdf
python3 scripts/render_figure4a_d6b.py
```

The manuscript figure command writes uncompressed 600-d.p.i. TIFF files to `outputs/main/`. The downstream SCZ figure command above writes the submission-format PDF with vector text and linework to `outputs/main/`; omit `--pdf` to produce the default uncompressed 600-d.p.i. TIFF. The output directories are separate from the submission package.

## Analysis path

1. Public GWAS files and phenotype representations are identified in the provenance tables.
2. `scripts/upstream/core/` standardizes source files and prepares LDSC-compatible inputs.
3. `scripts/run_ldsc_h2.sh` and `scripts/run_ldsc_rg.sh` run the LDSC h2 and rg commands against an external LDSC installation and LD reference.
4. `scripts/upstream/ascertainment_pilot/` computes PHBC, block-jackknife, leave-one-out, neurological and SCZ sensitivity results.
5. `scripts/upstream/p26_broad_scz/` computes the broad FinnGen SCZ representation sensitivity analysis.
6. `scripts/upstream/ascertainment_pilot/prepare_supergnova_inputs.py` prepares local-covariance inputs. `scripts/upstream/supergnova/run_supergnova.sh` launches an external SUPERGNOVA installation pinned in `environment/external_software.tsv`.
7. The final primary SCZ downstream Genomic SEM analysis is implemented in `scripts/upstream/ascertainment_pilot/d5i_joint_fullsample.R`. It uses the joint LDSC S/V input, fits the shared reduced compulsive, internalizing and substance-use domain structure with PGC and FinnGen SCZ as external traits, and writes the full-sample fit, contrast and boundary-sensitivity audit to `results/d5i_joint_fullsample/` when the corresponding remote analysis inputs are available. The targeted revision check removed ANX from COMP while retaining ANX in INT and is archived in `results/d5i_targeted_alt_20260905/`.
8. The D6A Stratified Genomic SEM analysis is implemented in `scripts/upstream/ascertainment_pilot/d6a_baseline_stratified_gsem.R`. It uses the official BaselineLD v2.2 annotations, the D5i-matched SCZ/domain input and `enrich()` with the LDSC 200-block sampling covariance. The referee-accessible long table is Supplementary Table 27.
9. The D6B published PI/brain-cell Stratified Genomic SEM analysis is archived under `results/d6b_published_pi_braincell_stratified_gsem_20260905/`. Its full output, custom subset and provenance/QC summary are Supplementary Tables 28–30. The covariance-aware direct PGC–FinnGen annotation contrasts are archived under `results/d6b_direct_contrast_20260905/` and reported in Supplementary Table 31. `scripts/upstream/ascertainment_pilot/d6b_direct_annotation_contrast.R` matches the six observed parameter keys, rescales the captured joint sandwich covariance and reconstructs the 492 FinnGen-minus-PGC contrasts. The targeted exploratory raw covariance pilot is reconstructed by `scripts/upstream/ascertainment_pilot/d6b_raw_covariance_contrast.R` for the fixed 29 PI/brain-cell annotations and three domains, with separate 29-test corrections and model/smoothing QC fields in Supplementary Table 32. `scripts/build_d6b_submission_tables.py` materializes the provenance/QC table, `scripts/render_figure4_d6b.py` renders the complete main PI/brain-cell Figure 4 and `scripts/render_figure4a_d6b.py` renders the standalone Figure 4a panel.
10. The resulting tables under `data/derived/` and `data/supplementary_tables/` are the referee-accessible outputs, including the primary SCZ Genomic SEM analysis with full parameter reporting in Supplementary Table 17, the six-trait SUPERGNOVA summary in Supplementary Table 21, the LDSC diagnostic and overlap audit in Supplementary Table 22, the SCZ downstream joint genetic R2 and Shapley outputs in Supplementary Tables 23–26, the BaselineLD Stratified Genomic SEM enrichment output in Supplementary Table 27, the PI/brain-cell Stratified Genomic SEM outputs in Supplementary Tables 28–30, the enrichment-ratio direct contrasts in Supplementary Table 31 and the targeted raw covariance contrasts in Supplementary Table 32.
10. `scripts/upstream/ascertainment_pilot/p32_observed_shapley_r2_decomposition.R` computes the observed covariance-based joint genetic R2, exact Shapley contributions, paired delete-one-block uncertainty and the profile-level sensitivity tests. With the same ten-trait panel, the joint genetic R2 equals the uncorrected PHBC quadratic form before bias correction. The full matrices and block caches are external analysis inputs. To rerun against an analysis root containing the recorded `results/p11_scz_neff_audit/`, `results/p25/` and `results/p26_broad_scz/` directories, use `ANALYSIS_ROOT=/path/to/phenotype_ascertainment_pilot Rscript scripts/upstream/ascertainment_pilot/p32_observed_shapley_r2_decomposition.R`; set `RESULTS_DIR` to redirect the output archive.
11. `scripts/render_manuscript_figures.py` combines the committed derived tables into Figures 1–3 and Supplementary Figures 1–3. `scripts/render_figure4_d6b.py` renders the main PI/brain-cell Figure 4 from Supplementary Table 29; pass `--pdf` for the submission-format PDF. `scripts/render_figure4a_d6b.py` renders the standalone Figure 4a panel as an uncompressed 600-d.p.i. TIFF and can emit hidden PDF/SVG QA files when output paths are supplied. `scripts/render_figure4_downstream.py` renders the covariance-based joint genetic R2 and Shapley figure for Supplementary Figure 4. The SCZ local-covariance renderer is kept in `scripts/render_figure3_diagnostics.py` and is called by the main renderer.

## Upstream environment

The upstream scripts require R 4.5.0, Python 3.9 for the pinned SUPERGNOVA environment, LDSC 1.0.1 at Git commit `2fdeeb3b44379408794154993dbd6101b8946b7e`, `pleioh2g` 0.1.3, GenomicSEM 0.0.5, `mvtnorm` 1.3-3 and the packages listed in `environment/`. The local LDSC audit rerun used Python 3.14.3 with a documented compatibility copy of that LDSC source. Reference LD resources are excluded from this repository and must be supplied through the environment variables described in `data/input/README.md`. The local submission snapshot is recorded in `environment/submission_snapshot.tsv`; the workspace has no Git metadata, so no remote tag was created.

The output mapping and input paths are listed in `analysis_manifest.tsv`.

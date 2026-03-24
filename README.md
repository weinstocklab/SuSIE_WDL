# UKB Fine-Mapping WDL

This repository contains a DNAnexus/WDL workflow for fine-mapping genome-wide significant REGENIE hits with `susie_rss`. The workflow consumes a single REGENIE summary-statistics parquet containing all phenotypes, identifies significant loci, computes LD from UK Biobank PLINK2 blocks, runs SuSiE per locus, and produces locus-level plots with gene tracks.

## Main Files

- [fine_mapping.WDL](/beegfs/labs/weinstocklab/projects/UKB/fine_mapping_WDL/fine_mapping.WDL): primary workflow
- [dx_compile.sh](/beegfs/labs/weinstocklab/projects/UKB/fine_mapping_WDL/dx_compile.sh): compile one or more WDLs with `dxCompiler`
- [create_job_submission.py](/beegfs/labs/weinstocklab/projects/UKB/fine_mapping_WDL/create_job_submission.py): render or submit `dx run` commands from YAML
- [fine_mapping_test_config.yaml](/beegfs/labs/weinstocklab/projects/UKB/fine_mapping_WDL/fine_mapping_test_config.yaml): example run config
- [job_submission.sh](/beegfs/labs/weinstocklab/projects/UKB/fine_mapping_WDL/job_submission.sh): saved submission command

## Required Inputs

The workflow `susie_finemapping` expects:

- `summary_stats_parquet`: concatenated REGENIE parquet with all phenotypes
- `step2_chunk_manifest`: manifest mapping loci to `pvar/psam/pgen` blocks
- `covariates_file`: GWAS covariate file used to define the analysis sample
- `plink2_binary`: zip containing the PLINK2 executable

Important behavior:

- Locus discovery applies `p_threshold` and `min_maf` (`0.001` by default).
- `covariates_file` is used only to build a PLINK keep file from `FID/IID`, `IID`, or `EID`.
- REGENIE chromosome `23` is normalized to `X`.
- LD computation handles loci with a single retained variant by emitting a `1x1` LD matrix.

## Outputs

Top-level outputs include:

- merged PIP table
- merged credible-set table
- per-locus `susie_rds`
- per-locus SuSiE logs
- per-locus locuszoom-style PDFs and PNGs
- per-locus plotting data tables
- loci manifest used for the scatter

## Typical Usage

Compile:

```bash
bash dx_compile.sh fine_mapping.WDL
```

Render a submission command:

```bash
python3 create_job_submission.py fine_mapping_test_config.yaml
```

Launch:

```bash
bash job_submission.sh
```

## Notes

- `compute_ld` uses UKB DRAGEN PLINK files and passes `--no-pheno` to PLINK2.
- Plotting is done in `plot_locus` and includes gene tracks based on `gwasplot`’s `gencode` annotations.
- If you change workflow inputs, recompile and then update both the YAML config and `job_submission.sh`.

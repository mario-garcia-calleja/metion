# Changelog

## 2.0.0

Complete rewrite of the Bachelor's thesis version (v1).

### Methods
- **Modification calling in all sequence contexts.** v1 used the `5mCG_5hmCG` model and
  `modkit pileup --cpg`, which restricts the output to CpG sites. Bacterial methyltransferases
  mostly act outside CpG (Dam: 6mA at `GATC`; Dcm: 5mC at `CCWGG`), so 6mA was effectively
  discarded and Dcm 5mC could not be detected. v2 uses the `6mA` and `4mC_5mC` all-context
  models and no CpG filter by default.
- **Statistical filtering.** v1 kept sites with `score > 0` from `modkit dmr pair`. v2 filters on
  minimum coverage, minimum change in methylated fraction and Benjamini–Hochberg FDR on the
  per-site MAP p-value. `--legacy` reproduces the v1 filter exactly.
- **Promoters.** Sites upstream of a region (`--flank`, strand-aware) are now reported.

### Engineering
- Single workflow (v1 needed two scripts and a hand-edited CSV of result paths).
  Samplesheet-driven: start from pod5, BAM or bedMethyl.
- Organism-agnostic: reference and regions BED are parameters.
- Failures stop the run (v1 silently ignored failed barcodes).
- Resource limits and `conda`/`mamba`/`micromamba`/`slurm`/`test` profiles.
- Analysis scripts rewritten with argparse; TSV outputs plus one Excel report; unit tests.
- Synthetic test dataset and GitHub Actions CI.

# Predicting Genetic Variants with Discordant Post-Imputation Allele Frequency

Supplementary code for the master's thesis of M. Filippov, ITMO University, 2023–2025.

Genotype imputation can produce variants whose allele frequencies differ from those observed in
sequencing data, which leads to false-positive associations. This project:

1. Runs QC on imputed 1000 Genomes genotypes, annotates them with Ensembl VEP and keeps protein-coding variants.
2. Clusters samples by ancestry and compares imputed and genotyped allele counts against
   ancestry-matched public exome controls from the [SCoRe](https://dnascore.net/) platform.
3. Flags variants with discordant allele frequencies (weighted linear regression, Fisher and χ² tests, BH-FDR).
4. Trains gradient-boosting classifiers (XGBoost, LightGBM, CatBoost, tuned with Optuna) to predict
   discordant variants from imputation quality, frequency, conservation and functional annotations,
   validated with leave-one-chromosome-out cross-validation.

## Repository structure

| File | Description |
|---|---|
| `QC_and_filtering.sh` | QC, VEP annotation, coding-variant filtering, typed/imputed split, per-cluster genotype counts |
| `SCoRe.R` | Ancestry PCA and clustering ([SVDFunctions](https://github.com/alexloboda/SVDFunctions)); writes the SCoRe input `1000G.yaml` |
| `GWAS.R` | Case–control association tests against SCoRe controls; writes `results_full_annotated.csv` |
| `exploration.Rmd` | Plots: MAF/R² distributions, Manhattan plot, p-value relationships |
| `training.ipynb` | Model selection, tuning and leave-one-chromosome-out evaluation |
| `tests/` | Synthetic data generator and an end-to-end smoke test |

## Requirements

Tested with bcftools 1.19, Ensembl VEP v113 (offline cache), R 4.4 and Python 3.10.

- R packages: `SVDFunctions` (`devtools::install_github("alexloboda/SVDFunctions")`), `dplyr`, `stringr`,
  plus `ggplot2`, `ggman`, `flextable`, `officer` for `exploration.Rmd`
- Python packages: see `requirements.txt`

## How to run

Run all commands from the repository folder.

1. Put the imputed VCF (Minimac-style `R2` and `TYPED` INFO fields) in the folder as `1000G.imputed.vcf.gz`.
2. Run QC and clustering:
   ```bash
   bash QC_and_filtering.sh 1000G.imputed.vcf.gz
   ```
   This writes `1000G.yaml`, `vars_typed.txt`, `vars_imputed.txt`, `variant_annotations.csv`
   and per-cluster case counts in `case_counts/`.
3. On the SCoRe platform, submit `1000G.yaml` with `vars_typed.txt` (default QC filters).
   Unpack the results into `controls_counts/` and rename each `counts_X.tsv` to `counts_typed_X.tsv`.
4. Repeat step 3 with `vars_imputed.txt`, renaming the files to `counts_imputed_X.tsv`
   (same `controls_counts/` folder).
5. Run the association analysis:
   ```bash
   Rscript GWAS.R
   ```
   The result is `results_full_annotated.csv`; explore it with `exploration.Rmd`.
6. Model training is in `training.ipynb`. No data is included in this repository. The notebook
   expects a `data.csv` with the columns listed in its first cell (the GWAS results joined with
   VEP plugin annotations: GERP, gnomAD, SIFT, LCR). Set `DATA_PATH` and `N_TRIALS` environment
   variables to override the input path and the number of Optuna trials.

## Testing

`tests/run_tests.sh` runs the whole pipeline end to end on synthetic data. It replaces the two
external dependencies (the VEP cache and the SCoRe web platform) with local mocks:

```bash
bash tests/run_tests.sh
```

The synthetic data is only for checking that the code runs. It does not reproduce thesis results.

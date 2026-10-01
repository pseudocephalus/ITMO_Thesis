#!/usr/bin/env bash
# End-to-end smoke test on synthetic data.
# Replaces the two external steps (Ensembl VEP, SCoRe web platform) with local mocks.
# Usage: bash tests/run_tests.sh
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
[[ -n "${KEEP:-}" ]] && echo "Working dir: $WORK" || trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
cp "$REPO"/QC_and_filtering.sh "$REPO"/SCoRe.R "$REPO"/GWAS.R .
mkdir -p bin && cp "$REPO/tests/mock_vep" bin/vep && chmod +x bin/vep
export PATH="$WORK/bin:$PATH"

echo "== 1. Mock imputed VCF"
Rscript -e 'd <- SVDFunctions::publicExomesDataset; write.table(data.frame(d$variants, d$mean, d$U[,1:3]), "score_variants.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)'
python3 "$REPO/tests/make_mock_data.py" vcf --score-variants score_variants.tsv

echo "== 2. QC, annotation, ancestry clustering"
bash QC_and_filtering.sh 1000G.imputed.vcf.gz
for f in vars_typed.txt vars_imputed.txt 1000G.yaml variant_annotations.csv 1000G_ids1.txt; do
    [[ -s "$f" ]] || { echo "FAIL: $f missing or empty"; exit 1; }
done
ls case_counts

echo "== 3. Mock SCoRe controls + association tests"
python3 "$REPO/tests/make_mock_data.py" controls
Rscript GWAS.R
[[ -s results_full_annotated.csv ]] || { echo "FAIL: no GWAS results"; exit 1; }
Rscript -e 'df <- read.csv("results_full_annotated.csv"); print(table(df$type, df$cluster)); cat("Significant (BH<0.05):", sum(df$p_regr_adj < 0.05), "of", nrow(df), "\n")'

echo "== 4. Exploration report"
if Rscript -e 'quit(status = !requireNamespace("ggman", quietly = TRUE))' 2>/dev/null; then
    cp "$REPO/exploration.Rmd" . && Rscript -e 'rmarkdown::render("exploration.Rmd", output_format = "html_document", quiet = TRUE)'
    echo "exploration.Rmd rendered OK"
else
    echo "skipped (R package ggman not installed)"
fi

echo "== 5. Model training notebook (mock features, 2 Optuna trials per model)"
python3 "$REPO/tests/make_mock_data.py" training --out "$WORK/data.csv"
cp "$REPO/training.ipynb" .
DATA_PATH="$WORK/data.csv" N_TRIALS=2 jupyter nbconvert --to notebook --execute \
    training.ipynb --output training_executed.ipynb \
    --ExecutePreprocessor.timeout=1800 >/dev/null 2>&1
echo "Notebook executed OK"

echo "ALL TESTS PASSED"

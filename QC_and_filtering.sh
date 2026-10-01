#!/usr/bin/env bash
# QC, annotation and filtering of imputed genotypes.
# Usage: bash QC_and_filtering.sh [imputed.vcf.gz]
set -euo pipefail

# Define variables
INPUT_VCF="${1:-1000G.imputed.vcf.gz}"
FILTERED_VCF="1000G.imputed.filtered.vcf.gz"
VEP_VCF="1000G.filtered.vep.vcf.gz"
EXOME_VCF="1000G.imputed.filtered.vep.exome.vcf.gz"
TYPED_VCF="1000G.imputed.filtered.vep.exome.typed.vcf.gz"
IMPUTED_VCF="1000G.imputed.filtered.vep.exome.imputed.vcf.gz"
EUR_SAMPLES="1kg_eur.txt"
CSQ_TERMS=(
    splice_acceptor_variant splice_donor_variant stop_gained
    frameshift_variant stop_lost start_lost inframe_insertion
    inframe_deletion missense_variant protein_altering_variant
    synonymous_variant
)

for tool in bcftools vep Rscript; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found in PATH" >&2; exit 1; }
done
[[ -f "$INPUT_VCF" ]] || { echo "ERROR: input VCF '$INPUT_VCF' not found" >&2; exit 1; }

# Step 1: Initial filtering (keep genotyped variants and well-imputed ones, R2 > 0.69)
bcftools view -i 'AF>0 && (R2>0.69 || TYPED==1)' "$INPUT_VCF" -Oz -o "$FILTERED_VCF"

# Step 2: VEP annotation
vep --offline --cache --pick --vcf --compress_output gzip --force_overwrite \
    -i "$FILTERED_VCF" -o "$VEP_VCF"

# Step 3: Create CSQ filter expression: CSQ~"term1" || CSQ~"term2" || ...
CSQ_FILTER=""
for term in "${CSQ_TERMS[@]}"; do
    CSQ_FILTER+="${CSQ_FILTER:+ || }CSQ~\"${term}\""
done

# Step 4: Filter for exome (protein-coding consequence) variants
bcftools view -i "$CSQ_FILTER" "$VEP_VCF" -Oz -o "$EXOME_VCF"

# Step 5: Split into typed/imputed
bcftools view -i 'TYPED=1' "$EXOME_VCF" -Oz -o "$TYPED_VCF"
bcftools view -e 'TYPED=1' "$EXOME_VCF" -Oz -o "$IMPUTED_VCF"

# Step 6: Extract variant lists (input for the SCoRe platform)
for type in typed imputed; do
    vcf_var="${type^^}_VCF"
    bcftools query -f '%CHROM:%POS %REF %ALT\n' "${!vcf_var}" > "vars_${type}.txt"
done

# Step 7: Download the list of European 1000 Genomes samples
[[ -f "$EUR_SAMPLES" ]] || wget -q http://dnascore.net/tutorial/1kg_eur.txt -O "$EUR_SAMPLES"

# Step 8: Ancestry clustering and SCoRe input (writes 1000G.yaml and 1000G_ids<N>.txt)
Rscript SCoRe.R "$VEP_VCF" "$EUR_SAMPLES"

# Step 9: Extract variant annotations
bcftools query -f '%ID,chr%CHROM,%POS,%REF,%ALT,%R2,%MAF,%CSQ\n' "$EXOME_VCF" > variant_annotations.csv

# Steps 10-11: Split into ancestry clusters and count genotypes per variant
mkdir -p case_counts
for ids in 1000G_ids*.txt; do
    cluster="${ids#1000G_ids}"
    cluster="${cluster%.txt}"
    for type in imputed typed; do
        in_var="${type^^}_VCF"
        out_vcf="1000G.imputed.filtered.vep.exome.${type}.cluster${cluster}.vcf.gz"
        bcftools view -S "$ids" "${!in_var}" -Oz -o "$out_vcf"

        # Columns: CHROM POS ID REF ALT unknown hom_ref het hom_alt alt_carriers
        bcftools view -H "$out_vcf" | awk '{
            unknown=homref=het=homalt=0
            for (i=10; i<=NF; i++) {
                split($i, a, ":")
                split(a[1], GT, "[/|]")
                if (GT[1]=="." && GT[2]==".") unknown++
                else if (GT[1]==0 && GT[2]==0) homref++
                else if (GT[1]==GT[2]) homalt++
                else het++
            }
            print $1, $2, $3":"$4"-"$5, $4, $5, unknown, homref, het, homalt, het+homalt
        }' > "case_counts/counts_${type}_${cluster}"
    done
done

echo "Done. Next: submit 1000G.yaml with vars_typed.txt / vars_imputed.txt to SCoRe (see README)."

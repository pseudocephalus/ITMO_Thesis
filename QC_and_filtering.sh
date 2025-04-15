#!/bin/bash

# Define variables
INPUT_VCF="1000G.imputed.vcf.gz"
FILTERED_VCF="1000G.imputed.filtered.vcf.gz"
VEP_VCF="1000G.filtered.vep.vcf.gz"
EXOME_VCF="1000G.imputed.filtered.vep.exome.vcf.gz"
TYPED_VCF="1000G.imputed.filtered.vep.exome.typed.vcf.gz"
IMPUTED_VCF="1000G.imputed.filtered.vep.exome.imputed.vcf.gz"
CSQ_TERMS=(
    splice_acceptor_variant splice_donor_variant stop_gained 
    frameshift_variant stop_lost start_lost inframe_insertion 
    inframe_deletion missense_variant protein_altering_variant 
    synonymous_variant
)
CLUSTERS=(1 2 3)

# Step 1: Initial filtering
bcftools view -i 'AF>0 && (R2>0.69 || TYPED==1)' "$INPUT_VCF" -Oz -o "$FILTERED_VCF"

# Step 2: VEP annotation
vep --offline --cache --pick --vcf --compress_output gzip \
    -i "$FILTERED_VCF" -o "$VEP_VCF"

# Step 3: Create CSQ filter expression
CSQ_FILTER=$(IFS="|"; echo "CSQ ~ \"${CSQ_TERMS[*]/#/\" || CSQ ~ \"}\"")
CSQ_FILTER="${CSQ_FILTER%||*}"  # Remove trailing ||

# Step 4: Filter for exome variants
bcftools view "$VEP_VCF" -i "$CSQ_FILTER" -Oz -o "$EXOME_VCF"

# Step 5: Split into typed/imputed
bcftools view -i 'TYPED=1' "$EXOME_VCF" -Oz -o "$TYPED_VCF"
bcftools view -e 'TYPED=1' "$EXOME_VCF" -Oz -o "$IMPUTED_VCF"

# Step 6: Extract variant lists
for type in typed imputed; do
    bcftools query -f '%CHROM:%POS %REF %ALT\n' "${type^^}_VCF" > "vars_${type}.txt"
done

# Step 7: Download population data
wget -q http://dnascore.net/tutorial/1kg_eur.txt

# Step 8: Run SCoRe.R
SCoRe.R

# Step 9: Extract variant annotations
bcftools query "$EXOME_VCF" -f '%ID,chr%CHROM,%POS,%REF,%ALT,%R2,%MAF,%CSQ' > variant_annotations.csv

# Step 10: Split into clusters
for cluster in "${CLUSTERS[@]}"; do
    for type in imputed typed; do
        in_vcf="${type^^}_VCF"
        out_vcf="1000G.imputed.filtered.vep.exome.${type}.cluster${cluster}.vcf.gz"
        bcftools view -S "1000G_ids${cluster}.txt" "${!in_vcf}" -Oz -o "$out_vcf"
    done
done

# Step 11: Generate counts
mkdir -p case_counts
for type in imputed typed; do
    for cluster in "${CLUSTERS[@]}"; do
        vcf="1000G.imputed.filtered.vep.exome.${type}.cluster${cluster}.vcf.gz"
        output="case_counts/counts_${type}_${cluster}"
        
        zgrep -v "^#" "$vcf" | awk '{
            unknown=homref=het=homalt=0
            for(i=10; i<=NF; i++) {
                split($i, a, ":")
                split(a[1], GT, "[/|]")
                if(GT[1]=="." && GT[2]==".") unknown++
                else if(GT[1]==0 && GT[2]==0) homref++
                else if(GT[1]==GT[2]) homalt++
                else het++
            }
            print $1, $2, $3":"$4"-"$5, $4, $5, unknown, homref, het, homalt, het+homalt
        }' > "$output"
    done
done
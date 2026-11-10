configfile: "config.yaml"

wildcard_constraints:
    type="typed|imputed",
    cluster=r"\d+"

CSQ_FILTER = " || ".join(f'CSQ~"{term}"' for term in config["csq_terms"])


rule all:
    input:
        "1000G.yaml",
        expand("vars_{type}.txt", type=config["types"]),
        "results_full_annotated.csv"


rule download_eur_samples:
    output:
        "1kg_eur.txt"
    params:
        url=config["eur_samples_url"]
    shell:
        "wget -q {params.url} -O {output}"


rule filter_qual:
    input:
        config["input_vcf"]
    output:
        "1000G.imputed.filtered.vcf.gz"
    params:
        r2=config["r2_threshold"]
    shell:
        "bcftools view -i 'AF>0 && (R2>{params.r2} || TYPED==1)' {input} -Oz -o {output}"


rule annotate:
    input:
        "1000G.imputed.filtered.vcf.gz"
    output:
        "1000G.filtered.vep.vcf.gz"
    shell:
        "vep --offline --cache --pick --vcf --compress_output gzip --force_overwrite "
        "-i {input} -o {output}"


rule filter_exome:
    input:
        "1000G.filtered.vep.vcf.gz"
    output:
        "1000G.imputed.filtered.vep.exome.vcf.gz"
    params:
        csq=CSQ_FILTER
    shell:
        "bcftools view -i '{params.csq}' {input} -Oz -o {output}"


rule split:
    input:
        "1000G.imputed.filtered.vep.exome.vcf.gz"
    output:
        typed="1000G.imputed.filtered.vep.exome.typed.vcf.gz",
        imputed="1000G.imputed.filtered.vep.exome.imputed.vcf.gz"
    shell:
        """
        bcftools view -i 'TYPED=1' {input} -Oz -o {output.typed}
        bcftools view -e 'TYPED=1' {input} -Oz -o {output.imputed}
        """


rule variant_list:
    input:
        "1000G.imputed.filtered.vep.exome.{type}.vcf.gz"
    output:
        "vars_{type}.txt"
    shell:
        r"bcftools query -f '%CHROM:%POS %REF %ALT\n' {input} > {output}"


rule score_clustering:
    input:
        vcf="1000G.filtered.vep.vcf.gz",
        samples="1kg_eur.txt",
        script="scripts/SCoRe.R"
    output:
        "1000G.yaml",
        expand("1000G_ids{cluster}.txt", cluster=config["clusters"])
    shell:
        "Rscript {input.script} {input.vcf} {input.samples}"


rule get_annotations:
    input:
        "1000G.imputed.filtered.vep.exome.vcf.gz"
    output:
        "variant_annotations.csv"
    shell:
        r"bcftools query -f '%ID,chr%CHROM,%POS,%REF,%ALT,%R2,%MAF,%CSQ\n' {input} > {output}"


rule cluster_subset:
    input:
        vcf="1000G.imputed.filtered.vep.exome.{type}.vcf.gz",
        ids="1000G_ids{cluster}.txt"
    output:
        "1000G.imputed.filtered.vep.exome.{type}.cluster{cluster}.vcf.gz"
    shell:
        "bcftools view -S {input.ids} {input.vcf} -Oz -o {output}"


rule count_genotypes:
    input:
        vcf="1000G.imputed.filtered.vep.exome.{type}.cluster{cluster}.vcf.gz",
        script="scripts/count_genotypes.awk"
    output:
        "case_counts/counts_{type}_{cluster}"
    shell:
        "bcftools view -H {input.vcf} | awk -f {input.script} > {output}"


rule gwas:
    input:
        cases=expand("case_counts/counts_{type}_{cluster}",
                     type=config["types"], cluster=config["clusters"]),
        controls=expand("controls_counts/counts_{type}_{cluster}.tsv",
                        type=config["types"], cluster=config["clusters"]),
        annotations="variant_annotations.csv",
        script="scripts/GWAS.R"
    output:
        "results_full_annotated.csv"
    shell:
        "Rscript {input.script}"
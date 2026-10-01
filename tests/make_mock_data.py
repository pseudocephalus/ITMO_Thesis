#!/usr/bin/env python3
"""Generate small synthetic inputs to smoke-test the pipeline (no real genotypes).

Subcommands:
  vcf       mock imputed VCF + sample list  (needs score_variants.tsv: id<TAB>mean)
  controls  mock SCoRe control counts       (run after QC_and_filtering.sh)
  training  mock feature table for training.ipynb
"""
import argparse
import glob
import os
import subprocess

import numpy as np
import pandas as pd

RNG = np.random.default_rng(666)
N_PER_POP = 100
N_POPS = 3


def make_vcf(args):
    # SCoRe variant ids are "chr1:123<TAB>REF<TAB>ALT", so each id spans three columns
    ref = pd.read_csv(args.score_variants, sep="\t", header=None,
                      names=["chr_pos", "ref", "alt", "mean", "u1", "u2", "u3"])
    ref["chrom"] = ref["chr_pos"].str.split(":").str[0].str.replace("chr", "", regex=False)
    ref["pos"] = ref["chr_pos"].str.split(":").str[1].astype(int)
    ref = ref[ref["chrom"].str.isdigit()].copy()
    ref["chrom_n"] = ref["chrom"].astype(int)
    ref = ref.sort_values(["chrom_n", "pos"]).reset_index(drop=True)
    n = len(ref)

    # three sub-populations shifted along the reference principal axes
    shifts = [ref["u1"] * 15, ref["u2"] * 15, -ref["u1"] * 15]
    samples, genos = [], []
    for k in range(N_POPS):
        af = np.clip((ref["mean"] + shifts[k]) / 2, 0.01, 0.99).to_numpy()
        g = RNG.binomial(1, af[:, None], (n, N_PER_POP)) + RNG.binomial(1, af[:, None], (n, N_PER_POP))
        genos.append(g)
        samples += [f"HG{k}{i:04d}" for i in range(N_PER_POP)]
    g = np.hstack(genos)

    typed = RNG.random(n) < 0.3
    r2 = np.where(typed, 1.0, RNG.uniform(0.5, 1.0, n))
    af_all = g.mean(axis=1) / 2
    gt = np.array(["0|0", "0|1", "1|1"])

    with open("mock.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in sorted(ref["chrom_n"].unique()):
            f.write(f"##contig=<ID={c}>\n")
        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Alt allele frequency">\n')
        f.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">\n')
        f.write('##INFO=<ID=R2,Number=1,Type=Float,Description="Imputation R2">\n')
        f.write('##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Genotyped variant">\n')
        f.write('##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Imputed variant">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for i, r in ref.iterrows():
            af = af_all[i]
            flag = "TYPED" if typed[i] else "IMPUTED"
            info = f"AF={af:.4f};MAF={min(af, 1 - af):.4f};R2={r2[i]:.3f};{flag}"
            f.write(f"{r.chrom}\t{r.pos}\t{r.chrom}:{r.pos}:{r.ref}:{r.alt}\t{r.ref}\t{r.alt}\t.\tPASS\t{info}\tGT\t"
                    + "\t".join(gt[g[i]]) + "\n")
    subprocess.run(["bcftools", "view", "mock.vcf", "-Oz", "-o", args.out], check=True)
    os.remove("mock.vcf")
    with open("1kg_eur.txt", "w") as f:
        f.write("\n".join(samples) + "\n")


def make_controls(args):
    """Control genotype counts in the SCoRe output layout; ~10% of imputed variants get a shifted AF."""
    os.makedirs("controls_counts", exist_ok=True)
    cols = ["chr", "pos", "X", "ref", "alt", "unknown", "hom_ref", "het", "hom_alt", "alt_sum"]
    for path in sorted(glob.glob("case_counts/counts_*_*")):
        vtype, cl = os.path.basename(path).split("_")[1:]
        cases = pd.read_csv(path, sep=" ", header=None, names=cols)
        n_case = cases[["hom_ref", "het", "hom_alt"]].sum(axis=1)
        af = (cases["het"] + 2 * cases["hom_alt"]) / (2 * n_case)
        if vtype == "imputed":
            discordant = RNG.random(len(af)) < 0.1
            af = np.where(discordant, np.clip(af + 0.15, 0, 0.95), af)
        n = 500
        hom_alt = RNG.binomial(n, af ** 2)
        het = RNG.binomial(n - hom_alt, np.clip(2 * af * (1 - af) / np.clip(1 - af ** 2, 1e-9, None), 0, 1))
        hom_ref = n - hom_alt - het
        with open(f"controls_counts/counts_{vtype}_{cl}.tsv", "w") as f:
            f.write("# mock SCoRe control counts\n# cluster " + cl + "\n# n_controls " + str(n) + "\n")
            f.write('""\t"hom_ref"\t"het"\t"hom_alt"\n')
            for c, p, r, a, h0, h1, h2 in zip(cases.chr, cases.pos, cases.ref, cases.alt, hom_ref, het, hom_alt):
                f.write(f'"chr{c}:{p}\t{r}\t{a}"\t{h0}\t{h1}\t{h2}\n')


def make_training(args):
    n = 3000
    df = pd.DataFrame({
        "chr_pos_gt": [f"chr{c}:{p}:A G" for c, p in zip(RNG.integers(1, 23, n), RNG.integers(1e5, 2e8, n))],
        "MAF": RNG.uniform(0.01, 0.5, n),
        "R2": RNG.uniform(0.7, 1.0, n),
        "het.x": RNG.integers(0, 100, n),
        "hom_alt.x": RNG.integers(0, 50, n),
        "hom_ref.x": RNG.integers(50, 300, n),
        "platform": RNG.choice(["array_A", "array_B"], n),
        "cluster": RNG.integers(1, 4, n),
        "GERP_RS": RNG.normal(2, 3, n),
        "GERP_NR": RNG.uniform(3, 6, n),
        "GERP_RS_rankscore": RNG.uniform(0, 1, n),
        "gnomADe_AF": RNG.uniform(0, 0.5, n),
        "gnomADg_AF": RNG.uniform(0, 0.5, n),
        "CSQ": RNG.choice(["missense_variant", "synonymous_variant", "stop_gained"], n),
        "BIOTYPE": RNG.choice(["protein_coding", "nonsense_mediated_decay"], n, p=[0.9, 0.1]),
        "SIFT": RNG.uniform(0, 1, n),
        "IMPACT": RNG.choice(["HIGH", "MODERATE", "LOW"], n),
        "LCR": RNG.random(n) < 0.05,
        "type": RNG.choice(["Genotyped", "Imputed"], n),
    })
    logit = -3 + 4 * (1 - df["R2"]) * 5 + 2 * df["LCR"] - 3 * df["MAF"]
    df["p_flag"] = np.where(RNG.random(n) < 1 / (1 + np.exp(-logit)), "Discordant", "Concordant")
    df.to_csv(args.out, index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("vcf"); p.add_argument("--score-variants", required=True); p.add_argument("--out", default="1000G.imputed.vcf.gz")
    sub.add_parser("controls")
    p = sub.add_parser("training"); p.add_argument("--out", default="data.csv")
    a = ap.parse_args()
    {"vcf": make_vcf, "controls": make_controls, "training": make_training}[a.cmd](a)

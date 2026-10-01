# Ancestry clustering of case samples and preparation of the SCoRe input file.
# Detailed description of each step can be found in https://dnascore.net/ (tutorial)
# Usage: Rscript SCoRe.R [annotated.vcf.gz] [1kg_eur.txt]

library(SVDFunctions)

args <- commandArgs(trailingOnly = TRUE)
VCFname_1000G <- if (length(args) >= 1) args[1] else "1000G.filtered.vep.vcf.gz"
samples_file  <- if (length(args) >= 2) args[2] else "1kg_eur.txt"

variants <- SVDFunctions::publicExomesDataset$variants
samples <- scan(samples_file, what = character())

gmatrix_1000G <- genotypeMatrixVCF(vcf = VCFname_1000G,
                                   DP = 0,
                                   GQ = 0,
                                   variants = variants,
                                   samples = samples,
                                   predictMissing = TRUE,
                                   verbose = TRUE)

# Drop samples/variants with call rate < 95%
gmatrix_1000G_filter <- filterGmatrix(gmatrix = gmatrix_1000G$genotype,
                                      imputationResults = gmatrix_1000G$predicted)

casePCA_1000G <- gmatrixPCA(gmatrix_1000G_filter$genotype, components = 3,
                            referenceMean = publicExomesDataset$mean,
                            SVDReference = publicExomesDataset$U)

caseCl_1000G <- estimateCaseClusters(PCA = casePCA_1000G,
                                     plotBIC = TRUE,
                                     plotDendrogram = TRUE,
                                     clusters = 10,
                                     minClusters = 3)

n_clusters <- length(caseCl_1000G$classes)
for (i in seq_len(n_clusters)) {
  caseCl_1000G[i] <- paste0("EUR", i)
}

cases_kept_1000G <- prepareInstance(gmatrix = gmatrix_1000G_filter$genotype,
                                    imputationResults = gmatrix_1000G_filter$predicted,
                                    controlsU = publicExomesDataset$U,
                                    meanControl = publicExomesDataset$mean,
                                    outputFileName = "1000G.yaml",
                                    title = "1000G",
                                    clusters = caseCl_1000G,
                                    keptSamplesFile = NULL)

# Leaf clusters have ids 1..n_clusters; save their sample lists
for (i in seq_len(n_clusters)) {
  writeLines(cases_kept_1000G[[i]], paste0("1000G_ids", i, ".txt"))
}
message("Found ", n_clusters, " clusters; wrote 1000G.yaml and 1000G_ids1..", n_clusters, ".txt")

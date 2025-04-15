# Detailed description of each step can be found in https://dnascore.net/ (tutorial)

library(SVDFunctions)

VCFname_1000G <- "1000G.filtered.vep.FC.vcf.gz"
variants <- SVDFunctions::publicExomesDataset$variants
samples <- scan("1kg_eur.txt", what = character())
gmatrix_1000G <- genotypeMatrixVCF(vcf = VCFname_1000G,
                                      DP = 0,
                                      GQ = 0,
                                      variants = variants,
                                      samples = samples,
                                      predictMissing = TRUE,
                                      verbose = TRUE)

casePCA_1000G <- gmatrixPCA(gmatrix_1000G_filter$gmatrix, components = 3,
                               referenceMean = publicExomesDataset$mean, 
                               SVDReference = publicExomesDataset$U)

caseCl_1000G <- estimateCaseClusters(PCA = casePCA_1000G,
                                        plotBIC = TRUE,
                                        plotDendrogram = TRUE,
                                        clusters = 10, 
                                        minClusters = 3)

caseCl_1000G[1] <- "EUR1"
caseCl_1000G[2] <- "EUR2"
caseCl_1000G[3] <- "EUR3"

cases_kept_1000G <-   prepareInstance(gmatrix = gmatrix_1000G$gmatrix,
                                         imputationResults = gmatrix_1000G$imputationResults,
                                         controlsU = publicExomesDataset$U, 
                                         meanControl = publicExomesDataset$mean, 
                                         outputFileName = "1000G.yaml",
                                         title = "1000G", 
                                         clusters = caseCl_1000G)

write(cases_kept_1000G[[1]],'1000G_ids1.txt')
write(cases_kept_1000G[[2]],'1000G_ids2.txt')
write(cases_kept_1000G[[3]],'1000G_ids3.txt')

